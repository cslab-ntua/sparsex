#include <cstdlib>
#include <iostream>
#include <inttypes.h>
#include <fstream>
#include <string>
#include <vector>
#include <time.h>
#include "../spm.h"

#define MSEC_PER_SEC 1000L
#define USEC_PER_SEC (MSEC_PER_SEC*MSEC_PER_SEC)
#define NSEC_PER_SEC (MSEC_PER_SEC*MSEC_PER_SEC*MSEC_PER_SEC)

using namespace std;
using namespace csx;


static inline bool elem_cmp_less(const SpmCooElem &e0,
                                 const SpmCooElem &e1)
{
	int ret;
	ret = CooCmp(static_cast<CooElem>(e0), static_cast<CooElem>(e1));
	return (ret < 0);
}

void print_elapsed_time(struct timespec *ts_start,
                        struct timespec *ts_end)
{
    long nsec_start = (long) ts_start->tv_sec*NSEC_PER_SEC + ts_start->tv_nsec;
    long nsec_end = (long) ts_end->tv_sec*NSEC_PER_SEC + ts_end->tv_nsec;
    long nsec_elapsed = nsec_end - nsec_start;

    printf("Elapsed time: %ld.%ld\n", nsec_elapsed / NSEC_PER_SEC,
           nsec_elapsed % NSEC_PER_SEC);
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        cerr << argv[0] << ": too few arguments" << endl;
        cerr << "Usage: " << argv[0] << " mmf_file" << endl;
        exit(1);
    }

    ifstream ifs(argv[1], ifstream::in);
    if (!ifs) {
        perror("error while opening mmf file");
        exit(1);
    }

    // Read mmf header
    uint64_t nr_rows, nr_cols, nr_nzeros;
    ifs >> nr_rows >> nr_cols >> nr_nzeros;

    vector<SpmCooElem> elems;
    elems.reserve(nr_nzeros);
    uint64_t elems_cnt = 0;
    while (ifs.good()) {
        SpmCooElem ce;
        if (elems_cnt == nr_nzeros)
            break;
        ifs >> ce.y >> ce.x >> ce.val >> skipws;
        ++elems_cnt;
        elems.push_back(ce);
    }
    assert(elems.size() == nr_nzeros);
    ifs.close();

    // Collect split points
    uint64_t w_size;
    const char *window_size = getenv("WINDOW_SIZE");
    if (!window_size)
        w_size = nr_rows;
    else
        w_size = atol(window_size);

    vector<uint64_t> split_points;
    vector<SpmCooElem>::const_iterator iter;
    uint64_t row, row_prev;
    row = row_prev = 0;
    uint64_t row_cnt = 0;
    split_points.push_back(0);
    for (iter = elems.begin(), elems_cnt = 0;
         iter != elems.end(); ++iter, ++elems_cnt) {
        SpmCooElem ce = *iter;
        row = ce.y;
        if (row == row_prev)
            continue;
        
        // New row
        ++row_cnt;
        if (row_cnt > w_size) {
            // New window
            split_points.push_back(elems_cnt);
            row_cnt = 0;
        }

        row_prev = row;
    }
    assert(elems_cnt == nr_nzeros);
    split_points.push_back(nr_nzeros);

    cout << "# splits = " << split_points.size() - 1 << endl;
    
    // Shuffle elements randomly
    random_shuffle(elems.begin(), elems.end());

    // Time standard sort
    struct timespec ts_start, ts_end;
    cout << "Sorting ... " << endl;
    clock_gettime(CLOCK_REALTIME, &ts_start);
    sort(elems.begin(), elems.end(), elem_cmp_less);
    clock_gettime(CLOCK_REALTIME, &ts_end);
    cout << "Full sort: ";
    print_elapsed_time(&ts_start, &ts_end);

    // Shuffle elements again
    random_shuffle(elems.begin(), elems.end());

    // Time partial sort
    vector<uint64_t>::const_iterator witer;
    cout << "Sorting partially ... " << endl;
    clock_gettime(CLOCK_REALTIME, &ts_start);
    for (witer = split_points.begin();
         witer != split_points.end() - 1; ++witer) {
        vector<SpmCooElem>::iterator e0, ee;
        e0 = elems.begin() + *witer;
        ee = elems.begin() + *(witer + 1);
        sort(e0, ee, elem_cmp_less);
    }
    clock_gettime(CLOCK_REALTIME, &ts_end);
    cout << "Partial sort: ";
    print_elapsed_time(&ts_start, &ts_end);
    
    return 0;
}
