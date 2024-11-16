//
// Created by hongbo on 16/11/24.
//

#ifndef LPCVT_RECON_UTILS_H
#define LPCVT_RECON_UTILS_H
#include <geogram/basic/numeric.h>
#include <geogram/basic/memory.h>
#include "CompareTriangles.h"
#include <geogram/basic/algorithm.h>
namespace NMCVT_PC{
using namespace GEO;
/**
     * \brief sine/cosine table.
     * \details We keep a small table of sines and cosines for
     *  speeding up things a little bit.
     *  Table entries are as follows:
     *  - sincos_table[i][0] = sin(2*M_PI*i/(sincos_nb-1))
     *  - sincos_table[i][1] = cos(2*M_PI*i/(sincos_nb-1))
 */
static double sincos_table[10][2] = {
    {0,1},
    {0.642788,0.766044},
    {0.984808,0.173648},
    {0.866025,-0.5},
    {0.34202,-0.939693},
    {-0.34202,-0.939693},
    {-0.866025,-0.5},
    {-0.984808,0.173648},
    {-0.642788,0.766044},
    {-2.44929e-16,1}
};

struct OrientNormal {
    /**
	 * \brief OrientNormal constructor.
	 * \param[in] v_in the index of a point
	 * \param[in] dot_in the dot product between the (unit)
	 *  normal vector at \p v_in and the normal vector at
	 *  the point that initiated propagation to \p v_in.
     */
    OrientNormal(
        index_t v_in, double dot_in
        ) : v(v_in), dot(dot_in) {
    }

    /**
	 * \brief Compares two OrientNormal objects
	 * \retval true if \p rhs should be processed before this
	 *  OrientObject.
	 * \retval false otherwise.
     */
    bool operator<(const OrientNormal& rhs) const {
        return (::fabs(dot) < ::fabs(rhs.dot));
    }
    index_t v;
    double dot;
};

static void co3ne_split_triangles_list(
    vector<index_t>& triangles,
    vector<index_t>& good_triangles,
    vector<index_t>& not_so_good_triangles
) {
    index_t nb_triangles = triangles.size()/3;

    // Step 1: normalize vertices order
    for(index_t i=0; i<triangles.size(); i+=3) {
        index_t* ptr = &triangles[i];
        std::sort(ptr, ptr+3);
    }

    // Step 2: sort the triangles in lexicographic order
    vector<index_t> t_sort(nb_triangles);
    for(index_t t=0; t<nb_triangles; ++t) {
        t_sort[t] = t ;
    }
    CompareTriangles compare_triangles(triangles);
    GEO::sort(t_sort.begin(), t_sort.end(), compare_triangles);


    // Step 3: select the triangles that appear exactly 3 times
    index_t if1 = 0;
    while(if1 < nb_triangles) {
        index_t if2 = if1 + 1;
        while(
            if2 < nb_triangles &&
            compare_triangles.is_same(t_sort[if1], t_sort[if2])
        ) {
            if2++;
        }

        index_t t = t_sort[if1];
        if(if2 - if1 == 3) {
            good_triangles.push_back(triangles[3*t]);
            good_triangles.push_back(triangles[3*t+1]);
            good_triangles.push_back(triangles[3*t+2]);
        } else if(if2 - if1 <= 2) {
            not_so_good_triangles.push_back(triangles[3*t]);
            not_so_good_triangles.push_back(triangles[3*t+1]);
            not_so_good_triangles.push_back(triangles[3*t+2]);
        }
        if1 = if2;
    }
}

}
#endif // LPCVT_RECON_UTILS_H
