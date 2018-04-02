#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <string>
#include <cstdint>
#include <unordered_map>

using namespace std;

unordered_map<string, int> nmf_bc_tbl({{"ONE_TO_ONE", 0}, {"SYM", 1}, {"WALL", 2}, {"FAR", 3}});

class NMF_BCRange {
public:
    uint32_t blk_seq;
    uint32_t face_seq;
    uint32_t pri_start_seq, pri_end_seq;
    uint32_t sec_start_seq, sec_end_seq;

    NMF_BCRange(uint32_t b = 0, uint32_t f = 0, uint32_t s1 = 0, uint32_t e1 = 0, uint32_t s2 = 0, uint32_t e2 = 0) :
            blk_seq(b), face_seq(f),
            pri_start_seq(s1), pri_end_seq(e1),
            sec_start_seq(s2), sec_end_seq(e2) {
        //Empty Initialization Body
    }
};

class NMF_Entry {
public:
    string bc;
    NMF_BCRange *rg1, *rg2;
    bool swap;

    NMF_Entry(string t) : bc(t), rg1(nullptr), rg2(nullptr), swap(false) {
        //Empty Initialization Body
    }
};

class NMF_BLKDim {
public:
    uint32_t i_dim, j_dim, k_dim;

    NMF_BLKDim(uint32_t in = 0, uint32_t jn = 0, uint32_t kn = 0) : i_dim(in), j_dim(jn), k_dim(kn) {
        //Empty Initialization Body
    }
};

class NMF {
public:
    uint32_t blk_num;
    vector<NMF_BLKDim> blk_dim;
    vector<NMF_Entry> mapping_entry;

    NMF(uint32_t n = 1) : blk_num(n), blk_dim(vector<NMF_BLKDim>(n)), mapping_entry(vector<NMF_Entry>()) {
        //Empty Initialization Body
    }
};

int main(int argc, char **argv) {
    if (argc != 3) {
        printf("Usage: gridglue map.nmf grid.fmt\n");
        return -1;
    }

    // Parse mapping file
    ifstream mfp;
    mfp.open(argv[1]);
    for (int i = 0; i < 4; i++)
        mfp.ignore(numeric_limits<streamsize>::max(), mfp.widen('\n'));


    uint32_t n = 0;
    mfp >> n;
    NMF bc_mp(n);
    for (uint32_t i = 0; i < n; i++) {
        uint32_t idx;
        mfp >> idx;
        mfp >> bc_mp.blk_dim[idx - 1].i_dim;
        mfp >> bc_mp.blk_dim[idx - 1].j_dim;
        mfp >> bc_mp.blk_dim[idx - 1].k_dim;
    }

    for (int i = 0; i < 4; i++)
        mfp.ignore(numeric_limits<streamsize>::max(), mfp.widen('\n'));

    string s;
    while (getline(mfp, s)) {
        cout << s << endl;
    }


    mfp.close();

    return 0;
}
