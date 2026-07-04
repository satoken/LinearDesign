#include <iomanip>
#include <random>
#include "beam_cky_parser.h"
#include "beam_cky_parser.cc"
#include "Utils/reader.h"
#include "Utils/common.h"
#include "Utils/codon.h"

// #ifndef CODON_TABLE
// #define CODON_TABLE "./codon_usage_freq_table_human.csv"
// #endif

#ifndef CODING_WHEEL
#define CODING_WHEEL "./coding_wheel.txt"
#endif

using namespace LinearDesign;

static string normalize_fixed_prefix(string fixed_prefix) {
    transform(fixed_prefix.begin(), fixed_prefix.end(), fixed_prefix.begin(), ::toupper);
    for (auto& nuc : fixed_prefix) {
        if (nuc == 'T')
            nuc = 'U';
        if (nuc != 'A' && nuc != 'C' && nuc != 'G' && nuc != 'U') {
            throw runtime_error("fixed RNA prefix must contain only A, C, G, U, or T");
        }
    }
    return fixed_prefix;
}

static pair<int, int> parse_avoid_pair_range(const string& range) {
    auto sep = range.find('-');
    if (sep == string::npos)
        sep = range.find(':');
    if (sep == string::npos)
        throw runtime_error("avoid pair range must be START-END, e.g. 10-50");

    int start = stoi(range.substr(0, sep));
    int end = stoi(range.substr(sep + 1));
    if (start < 1 || end < start)
        throw runtime_error("avoid pair range must be 1-based and satisfy START <= END");
    return make_pair(start - 1, end - 1);
}

template <typename IndexType>
static std::unordered_map<string, Lattice<IndexType>> perturb_lattice_weights(
        const std::unordered_map<string, Lattice<IndexType>>& base_lattices,
        std::mt19937_64& rng, const double gumbel_loc,
        const double gumbel_scale, const double perturb_strength) {
    auto perturbed_lattices = base_lattices;
    if (perturb_strength == 0.0)
        return perturbed_lattices;

    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    for (auto& aa_lattice : perturbed_lattices) {
        for (auto& node_edges : aa_lattice.second.right_edges) {
            for (auto& edge : node_edges.second) {
                double u = uniform(rng);
                while (u <= 0.0)
                    u = uniform(rng);
                const double gumbel = gumbel_loc - gumbel_scale * log(-log(u));
                std::get<2>(edge) += perturb_strength * gumbel;
            }
        }
    }
    return perturbed_lattices;
}

template <typename ScoreType, typename IndexType>
bool output_result(const DecoderResult<ScoreType, IndexType>& result, 
        const double duration, const double lambda, const bool is_verbose, 
        const Codon& codon, string& CODON_TABLE, const size_t coding_start = 0,
        const bool uses_pair_penalty = false) {

    stringstream ss;
    const string coding_sequence = result.sequence.substr(coding_start);
    if (is_verbose)
        ss << "Using lambda = " << (lambda / 100.) << "; Using codon frequency table = " << CODON_TABLE << endl;
    ss << "mRNA sequence:  " << result.sequence << endl;
    ss << "mRNA structure: " << result.structure << endl;
    if (uses_pair_penalty)
        ss << "mRNA pseudo-adjusted folding free energy: ";
    else
        ss << "mRNA folding free energy: ";
    ss << std::setprecision(2) << fixed << result.score << " kcal/mol; mRNA CAI: "
       << std::setprecision(3) << fixed << codon.calc_cai(coding_sequence) << endl;
    if (is_verbose)
        ss << "Runtime: " << duration << " seconds" << endl;
    cout << ss.str() << endl;

    return true;
}

void show_usage() {
    cerr << "echo SEQUENCE | ./lineardesign -l [LAMBDA]" << endl;
    cerr << "OR" << endl;
    cerr << "cat SEQ_FILE_OR_FASTA_FILE | ./lineardesign -l [LAMBDA]" << endl;
    cerr << "Optional: --fixedprefix RNA_PREFIX fixes a coding RNA prefix before the input amino-acid sequence" << endl;
    cerr << "Optional: --fixedutr RNA_PREFIX fixes a 5' UTR RNA prefix excluded from translation and CAI" << endl;
    cerr << "Optional: --avoidpairrange START-END --avoidpairpenalty KCAL discourages base pairs involving that 1-based range" << endl;
    cerr << "Optional: --samples N --gumbelloc LOC --gumbelscale SCALE --perturbstrength W --seed SEED runs perturb-and-map sampling" << endl;
}


int main(int argc, char** argv) {

    // default args
    double lambda = 0.0f;
    bool is_verbose = false;
    string CODON_TABLE = "./codon_usage_freq_table_human.csv";
    string fixed_prefix;
    string fixed_utr;
    string avoid_pair_range;
    double avoid_pair_penalty = 0.0;
    int avoid_pair_start = -1;
    int avoid_pair_end = -1;
    int sample_count = 1;
    double gumbel_loc = 0.0;
    double gumbel_scale = 1.0;
    double perturb_strength = 0.0;
    long long random_seed = -1;

    // parse args
    if (argc < 4 || argc == 7 || (argc > 8 && argc < 12) || argc > 13) {
        show_usage();
        return 1;
    }else{
        lambda = atof(argv[1]);
        is_verbose = atoi(argv[2]) == 1;
        if (string(argv[3]) != ""){
            CODON_TABLE = argv[3];
        }
        if (argc == 5) {
            try {
                fixed_prefix = normalize_fixed_prefix(argv[4]);
            } catch (const exception& e) {
                cerr << e.what() << endl;
                return 1;
            }
        }
        if (argc == 6) {
            try {
                fixed_prefix = normalize_fixed_prefix(argv[4]);
                fixed_utr = normalize_fixed_prefix(argv[5]);
            } catch (const exception& e) {
                cerr << e.what() << endl;
                return 1;
            }
        }
        if (argc >= 8) {
            try {
                fixed_prefix = normalize_fixed_prefix(argv[4]);
                fixed_utr = normalize_fixed_prefix(argv[5]);
                avoid_pair_range = argv[6];
                avoid_pair_penalty = atof(argv[7]);
                if (!avoid_pair_range.empty()) {
                    auto parsed_range = parse_avoid_pair_range(avoid_pair_range);
                    avoid_pair_start = parsed_range.first;
                    avoid_pair_end = parsed_range.second;
                    if (avoid_pair_penalty < 0.0)
                        throw runtime_error("avoid pair penalty must be non-negative");
                }
            } catch (const exception& e) {
                cerr << e.what() << endl;
                return 1;
            }
        }
        if (argc >= 12) {
            sample_count = atoi(argv[8]);
            gumbel_loc = atof(argv[9]);
            gumbel_scale = atof(argv[10]);
            perturb_strength = atof(argv[11]);
            if (argc >= 13)
                random_seed = atoll(argv[12]);
            if (sample_count < 1) {
                cerr << "samples must be at least 1" << endl;
                return 1;
            }
            if (gumbel_scale <= 0.0) {
                cerr << "gumbel scale must be positive" << endl;
                return 1;
            }
            if (perturb_strength < 0.0) {
                cerr << "perturb strength must be non-negative" << endl;
                return 1;
            }
            if (random_seed < -1) {
                cerr << "seed must be -1 or non-negative" << endl;
                return 1;
            }
        }
    } 
    lambda *= 100.;
    
    // load codon table and coding wheel
    Codon codon(CODON_TABLE);
    string fixed_prefix_aa;
    if (!fixed_prefix.empty()) {
        if (fixed_prefix.length() % 3 != 0) {
            cerr << "fixed RNA prefix length must be a multiple of 3 for coding-prefix design" << endl;
            return 1;
        }
        try {
            fixed_prefix_aa = codon.cvt_rna_seq_to_aa_seq(fixed_prefix);
        } catch (const exception& e) {
            cerr << "fixed RNA prefix is not valid for the codon table: " << e.what() << endl;
            return 1;
        }
        if (fixed_prefix_aa.find('*') != string::npos) {
            cerr << "fixed RNA prefix must not contain a stop codon" << endl;
            return 1;
        }
    }
    std::unordered_map<string, Lattice<IndexType>> aa_graphs_with_ln_weights;
    std::unordered_map<std::string, std::unordered_map<std::tuple<NodeType, NodeType>, std::tuple<double, NucType, NucType>, std::hash<std::tuple<NodeType, NodeType>>>> best_path_in_one_codon_unit;
    std::unordered_map<std::string, std::string> aa_best_path_in_a_whole_codon;
    prepare_codon_unit_lattice<IndexType>(CODING_WHEEL, codon, aa_graphs_with_ln_weights, best_path_in_one_codon_unit, aa_best_path_in_a_whole_codon, lambda);
    std::random_device random_device;
    std::mt19937_64 rng(random_seed >= 0
            ? static_cast<std::mt19937_64::result_type>(random_seed)
            : static_cast<std::mt19937_64::result_type>(random_device()));

    // main loop
    string aa_seq, aa_tri_seq;
    vector<string> aa_seq_list, aa_name_list;
    // load input
    for (string seq; getline(cin, seq);){
        if (seq.empty()) continue;
        if (seq[0] == '>'){
            aa_name_list.push_back(seq); // sequence name
            if (!aa_seq.empty())
                aa_seq_list.push_back(aa_seq);
            aa_seq.clear();
            continue;
        }else{
            rtrim(seq);
            aa_seq += seq;
        }
    }
    if (!aa_seq.empty())
        aa_seq_list.push_back(aa_seq);

    // start design
    for(int i = 0; i < aa_seq_list.size(); i++){
        if (aa_name_list.size() > i)
            cout << aa_name_list[i] << endl;
        auto& aa_seq = aa_seq_list[i];
        // convert to uppercase
        transform(aa_seq.begin(), aa_seq.end(), aa_seq.begin(), ::toupper);
        aa_tri_seq.clear();
        if (is_verbose)
            cout << "Input protein: " << aa_seq << endl;
        string design_aa_seq = fixed_prefix_aa + aa_seq;
        if (is_verbose && !fixed_prefix.empty())
            cout << "Fixed RNA prefix: " << fixed_prefix << "; translated prefix protein: " << fixed_prefix_aa << endl;
        if (is_verbose && !fixed_utr.empty())
            cout << "Fixed 5' UTR: " << fixed_utr << endl;
        if (!ReaderTraits<Fasta>::cvt_to_seq(design_aa_seq, aa_tri_seq)) 
            continue;

        auto protein = util::split(aa_tri_seq, ' ');
        for (int sample_index = 0; sample_index < sample_count; ++sample_index) {
            auto sample_graphs_with_ln_weights = perturb_lattice_weights<IndexType>(
                    aa_graphs_with_ln_weights, rng, gumbel_loc, gumbel_scale, perturb_strength);
            if (sample_count > 1)
                cout << "# sample " << (sample_index + 1) << endl;

            // init parser
            BeamCKYParser<ScoreType, IndexType> parser(lambda, is_verbose);
            if (!avoid_pair_range.empty()) {
                parser.set_pair_penalty(avoid_pair_start, avoid_pair_end,
                        static_cast<ScoreType>(avoid_pair_penalty * 100.0));
                if (is_verbose && sample_index == 0)
                    cout << "Avoiding base pairs involving nt " << (avoid_pair_start + 1)
                     << "-" << (avoid_pair_end + 1) << " with pseudo-energy "
                     << avoid_pair_penalty << " kcal/mol per pair" << endl;
            }
            if (is_verbose && sample_count > 1 && sample_index == 0) {
                cout << "Perturb-and-map samples = " << sample_count
                     << "; Gumbel(loc=" << gumbel_loc << ", scale=" << gumbel_scale
                     << "); strength = " << perturb_strength
                     << "; seed = " << random_seed << endl;
            }

            // parse
            auto system_start = chrono::system_clock::now();
            auto dfa = get_dfa<IndexType>(sample_graphs_with_ln_weights, util::split(aa_tri_seq, ' '), fixed_prefix, fixed_utr);
            auto result = parser.parse(dfa, codon, design_aa_seq, protein, aa_best_path_in_a_whole_codon, best_path_in_one_codon_unit, sample_graphs_with_ln_weights, fixed_utr.length());
            auto system_diff = chrono::system_clock::now() - system_start;
            auto system_duration = chrono::duration<double>(system_diff).count();  

            // output
            output_result(result, system_duration, lambda, is_verbose, codon, CODON_TABLE,
                    fixed_utr.length(), !avoid_pair_range.empty() && avoid_pair_penalty > 0.0);

#ifdef FINAL_CHECK
            if (!fixed_utr.empty() && result.sequence.substr(0, fixed_utr.length()) != fixed_utr) {
                std::cerr << "Fixed UTR Check Failed:" << std::endl;
                std::cerr << result.sequence.substr(0, fixed_utr.length()) << std::endl;
                std::cerr << fixed_utr << std::endl;
                assert(false);
            }
            if (!fixed_prefix.empty() && result.sequence.substr(fixed_utr.length(), fixed_prefix.length()) != fixed_prefix) {
                std::cerr << "Fixed Prefix Check Failed:" << std::endl;
                std::cerr << result.sequence.substr(fixed_utr.length(), fixed_prefix.length()) << std::endl;
                std::cerr << fixed_prefix << std::endl;
                assert(false);
            }
            auto coding_sequence = result.sequence.substr(fixed_utr.length());
            if (codon.cvt_rna_seq_to_aa_seq(coding_sequence) != design_aa_seq) {
                std::cerr << "Final Check Failed:" << std::endl;
                std::cerr << codon.cvt_rna_seq_to_aa_seq(coding_sequence) << std::endl;
                std::cerr << design_aa_seq << std::endl;
                assert(false);
            }
#endif
        }
    }
    return 0;
}
