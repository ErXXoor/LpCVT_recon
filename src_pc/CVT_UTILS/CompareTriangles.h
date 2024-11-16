//
// Created by hongbo on 16/11/24.
//

#ifndef LPCVT_RECON_COMPARETRIANGLES_H
#define LPCVT_RECON_COMPARETRIANGLES_H
#include <geogram/basic/numeric.h>
#include <geogram/basic/memory.h>
namespace NMCVT_PC{
using namespace GEO;
class CompareTriangles {
  public:
    /**
         * \brief Constructs a new CompareFacets.
         * \param[in] triangles a const reference to a vector
         *  of indices triplets
     */
    explicit CompareTriangles(const vector<index_t>& triangles) :
                                                                  triangles_(triangles) {
    }

    /**
         * \brief Tests the lexicographic order of two facets by their indices.
         * \param[in] f1 index of the first facet
         * \param[in] f2 index of the second facet
         * \return true if facet \p f1 is before facet \p f2 according to
         *  the lexicographic order of its vertices, false otherwise.
     */
    bool is_before(index_t f1, index_t f2) const {
        for(index_t c=0; c<3; c++) {
            index_t v1 = triangles_[3*f1+c];
            index_t v2 = triangles_[3*f2+c];
            if(v1 > v2) {
                return false;
            }
            if(v1 < v2) {
                return true;
            }
        }
        return false;
    }

    /**
         * \brief Tests whether two facets are identical.
         * \param[in] f1 index of the first facet
         * \param[in] f2 index of the second facet
         * \return true if facets \p f1 and \p f2 have the same
         *  vertices, false otherwise
     */
    bool is_same(index_t f1, index_t f2) const {
        for(index_t c=0; c<3; c++) {
            index_t v1 = triangles_[3*f1+c];
            index_t v2 = triangles_[3*f2+c];
            if(v1 != v2) {
                return false;
            }
        }
        return true;
    }

    /**
         * \brief Tests the lexicographic order of two facets by their indices.
         * \param[in] f1 index of the first facet
         * \param[in] f2 index of the second facet
         * \return true if facet \p f1 is before facet \p f2 according to
         *  the lexicographic order of its vertices, false otherwise.
     */
    bool operator() (index_t f1, index_t f2) const {
        return is_before(f1, f2);
    }

  private:
    const vector<index_t>& triangles_;
};

}
#endif // LPCVT_RECON_COMPARETRIANGLES_H
