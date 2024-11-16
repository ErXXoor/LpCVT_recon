//
// Created by hongbo on 16/11/24.
//

#ifndef LPCVT_RECON_CO3NEMANIFOLDEXTRACTION_H
#define LPCVT_RECON_CO3NEMANIFOLDEXTRACTION_H
#include <geogram/basic/numeric.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/basic/numeric.h>
#include <geogram/basic/command_line_args.h>
#include <stack>
namespace NMCVT_PC{
using namespace GEO;
class Co3NeManifoldExtraction {
  public:
    static const index_t NO_CORNER = index_t(-1);
    static const index_t NO_FACET  = index_t(-1);
    static const index_t NO_CNX    = index_t(-1);

    /**
         * \brief Initializes a new Co3NeManifoldExtraction with
         *  a list of triangles.
         * \param[in,out] target the target mesh. It needs to be already
         *  initialized with the vertices.
         * \param[in,out] good_triangles the good triangles reconstructed
         *  by Co3Ne. They are 'stealed' by the mesh (on exit, good_triangles
         *  is empty). If some non-manifold edges are detected, then all
         *  the triangles incident to any manifold edge are ignored.
     */
    Co3NeManifoldExtraction(
        Mesh& target,
        vector<index_t>& good_triangles
        ) : M_(target) {
        strict_ = false;
        if(strict_) {
            vector<index_t> first_triangle;
            for(index_t i=0; i<3; ++i) {
                first_triangle.push_back(*good_triangles.rbegin());
                good_triangles.pop_back();
            }

            M_.facets.assign_triangle_mesh(first_triangle, true);

            init_and_remove_non_manifold_edges();
            init_connected_components();
            add_triangles(good_triangles);
        } else {

            M_.facets.assign_triangle_mesh(good_triangles, true);

            init_and_remove_non_manifold_edges();
            init_connected_components();
        }
    }

    /**
         * \brief Tentatively adds triangle from the specified list.
         * \details Some geometric and topological properties are
         *  verified by connect_and_validate_triangle() before accepting
         *  the triangle.
         * \see connect_and_validate_triangle()
     */
    void add_triangles(const vector<index_t>& not_so_good_triangles) {
        bool pretty = false;

        index_t nb_triangles = not_so_good_triangles.size()/3;
        Logger::out("Co3ne") << "Tentatively add "
                             << nb_triangles << " triangles" << std::endl;
        vector<bool> t_is_classified(nb_triangles,false);
        bool changed = true;
        index_t max_iter = strict_ ? 5000 : 50;
        index_t iter = 0;
        bool first = true;
        while(changed && iter < max_iter) {
//            if(first) {
//                CmdLine::ui_clear_line();
//            } else {
//                first = false;
//            }
//            if(pretty) {
//                CmdLine::ui_message(
//                    "o-[Manifold Rec] Iteration:" + String::to_string(iter)
//                );
//            } else {
                Logger::out("Manifold Rec")
                    << "Iteration:" << iter << std::endl;
//            }
            changed = false;
            ++iter;
            for(index_t t=0; t<nb_triangles; ++t) {
                if(!t_is_classified[t]) {
                    index_t i = not_so_good_triangles[3*t];
                    index_t j = not_so_good_triangles[3*t+1];
                    index_t k = not_so_good_triangles[3*t+2];
                    index_t new_t = add_triangle(i,j,k);
                    bool classified = false;
                    if(connect_and_validate_triangle(new_t, classified)) {
                        changed = true;
                    } else {
                        rollback_triangle();
                    }
                    if(classified) {
                        t_is_classified[t] = true;
                    }
                }
            }
        }
//        if(pretty) {
//            CmdLine::ui_clear_line();
//            CmdLine::ui_message(
//                "o-[Manifold Rec] Iteration:" +
//                String::to_string(iter) + "\n"
//            );
//        } else {
//            Logger::out("Manifold Rec")
//                << "Iteration:" << iter << std::endl;
//        }
    }

  protected:

    /**
         * \brief Initializes the combinatorial data
         *  structures and deletes all facets incident
         *  to a non-manifold edge.
     */
    void init_and_remove_non_manifold_edges() {
        next_c_around_v_.assign(
            M_.facet_corners.nb(), index_t(NO_CORNER)
        );
        v2c_.assign(
            M_.vertices.nb(),index_t(NO_CORNER)
        );
        for(index_t t=0; t<M_.facets.nb(); ++t) {
            insert(t);
        }

        index_t nb_non_manifold = 0;
        vector<index_t> remove_t;
        for(index_t t=0; t<M_.facets.nb(); ++t) {
            if(!connect(t)) {
                remove_t.resize(M_.facets.nb(),0);
                remove_t[t] = 1;
                ++nb_non_manifold;
            }
        }

        mesh_reorient(M_, &remove_t);

        if(remove_t.size() == 0) {
            Logger::out("Co3Ne")
                << "All edges are manifold and well oriented"
                << std::endl;
        } else {
            index_t nb_remove_t = 0;
            for(index_t t=0; t<M_.facets.nb(); ++t) {
                if(remove_t[t] != 0) {
                    ++nb_remove_t;
                }
            }
            index_t nb_moebius = nb_remove_t - nb_non_manifold;
            Logger::out("Co3Ne")
                << "Removing " << nb_remove_t
                << " triangles ("
                << nb_non_manifold << " non_manifold, "
                << nb_moebius
                << " moebius)"
                << std::endl;
            M_.facets.delete_elements(remove_t,false);
        }

        // We need to re-compute next_c_around_v_ and v2c_
        // since all the indices changed in the mesh
        // (even if remove_t is empty, because mesh_reorient() may
        // have changed triangles orientation).
        next_c_around_v_.assign(M_.facet_corners.nb(), index_t(NO_CORNER));
        v2c_.assign(M_.vertices.nb(),index_t(NO_CORNER));
        for(index_t t=0; t<M_.facets.nb(); ++t) {
            insert(t);
        }
    }

    /**
         * \brief Tentatively connects a newly added triangle
         *  to the current mesh under construction. Accepted
         *  triangles satisfy the following criteria:
         *  - each new triangle should be either incident to at least
         *    two edges of existing triangles, or to one existing triangle
         *    and one isolated point.
         *  - the normals to the new triangle and its neighbor should
         *    not point to opposite directions.
         *  - inserting the new triangle should not generate 'by-excess'
         *    non-manifold vertices. A 'by-excess' non-manifold vertex
         *    has a closed loop of triangles in its neighbors plus
         *    additional triangles.
         *  - the orientation of the surface should be coherent (no Moebius
         *    strip).
         * \param[in] t index of the triangle
         * \param[out] classified true if the status of the triangle
         *  (accepted/rejected) could be completely determined,
         *  false if its status may still change during subsequent iterations
         * \retval true if all combinatorial and geometric tests succeeded
         * \retval false otherwise
     */
    bool connect_and_validate_triangle(index_t t, bool& classified) {
        index_t adj_c[3];
        classified = false;

        //   Combinatorial test (I): tests whether the three
        // candidate edges are manifold.
        if(!get_adjacent_corners(t,adj_c)) {
            classified = true;
            return false ;
        }

        //   Geometric test: tests whether the angles formed with
        // the candidate neighbors do not indicate degenerate sharp
        // creases.
        for(index_t i=0; i<3; ++i) {
            if(adj_c[i] != NO_CORNER) {
                index_t t2 = c2f(adj_c[i]);
                if(!triangles_normals_agree(t,t2)) {
                    classified = true;
                    return false;
                }
            }
        }

        int nb_neighbors =
            (adj_c[0] != NO_CORNER) +
            (adj_c[1] != NO_CORNER) +
            (adj_c[2] != NO_CORNER) ;

        // Combinatorial test (II)
        switch(nb_neighbors) {
            // If the candidate triangle is adjacent to no other
            // triangle, reject it
        case 0: {
            return false ;
        }
        // If the candidate triangle is adjacent to a single
        // triangle, reject it if the vertex opposite to
        // the common edge is not isolated.
        case 1: {
            // If not in strict mode, we reject the triangle.
            // Experimentally, it improves the result.
            if(!strict_) {
                return false;
            }
            index_t other_vertex=index_t(-1);
            for(index_t i=0; i<3; ++i) {
                if(adj_c[i] != NO_CORNER) {
                    other_vertex =
                        M_.facet_corners.vertex(
                            M_.facets.corners_begin(t) + ((i+2)%3)
                        );
                }
            }
            geo_debug_assert(other_vertex != index_t(-1));
            // Test whether other_vertex is isolated, reject
            // the triangle if other_vertex is NOT isolated.
            index_t nb_incident_T = nb_incident_triangles(other_vertex);
            geo_assert(nb_incident_T != 0); // There is at least THIS T.
            if(nb_incident_T > 1) {
                return false;
            }
        }
        }

        connect_adjacent_corners(t,adj_c);

        // Combinatorial test (III): test non-manifold vertices
        for(
            index_t c=M_.facets.corners_begin(t);
            c<M_.facets.corners_end(t); ++c
        ) {
            index_t v = M_.facet_corners.vertex(c);
            bool moebius=false;
            if(vertex_is_non_manifold_by_excess(v,moebius)) {
                classified = true;
                return false;
            }
            // It should not occur since we remove all Moebius configs
            // from the T3s and forbid Moebius configs when inserting
            // the T12s. However, some transient moebius configurations
            // due to triangle t may appear (since the Moebius test is
            // right after the non-manifold test).
            if(moebius) {
                Logger::warn("Co3Ne")
                    << "Encountered Moebius configuration" << std::endl;
                classified = true;
                return false;
            }
        }

        // Combinatorial test (IV): orientability
        if(!enforce_orientation_from_triangle(t)) {
            return false;
        }

        classified = true;
        return true;
    }


    /**
         * \brief Tentatively enforces mesh orientation starting from a
         *  given triangle.
         * \details The triangle \p t is rejected if it is incident to
         *  the same connected component with two different orientations.
         * \param[in] t index of the triangle to start mesh orientation from
         * \retval true if the mesh could be coherently oriented
         * \retval false otherwise
     */
    bool enforce_orientation_from_triangle(index_t t) {

        // Index of adjacent triangle
        // (or NO_FACET if no neighbor)
        index_t adj[3];

        // Index of adjacent connected component
        // (or NO_CNX if no neighbor)
        index_t adj_cnx[3];

        //   Orientation of adjacent triangle relative to
        // triangle t (or 0 if no neighbor)
        signed_index_t adj_ori[3];

        for(index_t i=0; i<3; ++i) {
            index_t c = M_.facets.corners_begin(t)+i;
            adj[i] = index_t(M_.facet_corners.adjacent_facet(c));
        }


        for(index_t i=0; i<3; ++i) {
            if(adj[i] == NO_FACET) {
                adj_ori[i] = 0;
                adj_cnx[i] = NO_CNX;
            } else {
                adj_ori[i] =
                    (triangles_have_same_orientation(t,adj[i])) ? 1 : -1;
                adj_cnx[i] = cnx_[adj[i]];
            }
        }

        //  If in the neighborhood the same connected component appears
        // with two opposite orientations, then connecting the triangle
        // would create a Moebius strip (the triangle is rejected)
        for(index_t i=0; i<3; ++i) {
            if(adj[i] != NO_FACET) {
                for(index_t j=i+1; j<3; ++j) {
                    if(
                        adj_cnx[j] == adj_cnx[i] &&
                        adj_ori[j] != adj_ori[i]
                    ) {
                        return false;
                    }
                }
            }
        }

        //  The triangle is accepted,
        // now reorient all the connected components and the
        // triangle coherently.

        // Find the largest component incident to t
        index_t largest_neigh_comp = NO_CNX;
        for(index_t i=0; i<3; ++i) {
            if(
                adj_cnx[i] != NO_CNX && (
                                            largest_neigh_comp == NO_CNX ||
                                            cnx_size_[adj_cnx[i]] >
                                                cnx_size_[adj_cnx[largest_neigh_comp]]
                                            )
            ) {

                largest_neigh_comp = i;
            }
        }
        geo_assert(largest_neigh_comp != NO_CNX);

        // Orient t like the largest incident component
        index_t comp = adj_cnx[largest_neigh_comp];

        cnx_.resize(std::max(t+1, cnx_.size()));
        cnx_[t] = comp;
        ++cnx_size_[comp];
        if(adj_ori[largest_neigh_comp] == -1) {
            flip_triangle(t);
            for(index_t i=0; i<3; ++i) {
                adj_ori[i] = -adj_ori[i];
            }
        }

        // Merge (and reorient if need be) all the other incident
        // components
        for(index_t i=0; i<3; ++i) {
            if(
                i != largest_neigh_comp &&
                adj[i] != NO_FACET && cnx_[adj[i]] != comp
            ) {
                merge_connected_component(
                    adj[i], comp, (adj_ori[i] == -1)
                );
            }
        }

        return true;
    }


    /**
         * \brief Adds a new triangle to the surface and to the
         *  combinatorial data structure.
         * \param[in] i first index of the triangle
         * \param[in] j second index of the triangle
         * \param[in] k third index of the triangle
     */
    index_t add_triangle(index_t i, index_t j, index_t k) {
        index_t result = M_.facets.create_triangle(i,j,k);
        next_c_around_v_.push_back(index_t(NO_CORNER));
        next_c_around_v_.push_back(index_t(NO_CORNER));
        next_c_around_v_.push_back(index_t(NO_CORNER));
        insert(result);
        return result;
    }

    /**
         * \brief Removes the latest triangle from both
         *  the mesh and the combinatorial data structure.
     */
    void rollback_triangle() {
        index_t t = M_.facets.nb()-1;
        remove(t);
        M_.facets.pop();
    }


    /**
         * \brief Inverts the orientation of a triangle.
         * \param[in] t the index of the triangle to be flipped.
     */
    void flip_triangle(index_t t) {

        // Remove t from the additional combinatorial data structure
        // (it is both simpler and more efficient to do that
        //  than updating it).
        remove(
            t,
            false // disconnect is set to false because
                  // we will re-insert t right after.
        );

        index_t c1 = M_.facets.corners_begin(t);
        index_t c2 = c1+1;
        index_t c3 = c2+1;
        index_t v1 = M_.facet_corners.vertex(c1);
        index_t f1 = M_.facet_corners.adjacent_facet(c1);
        index_t f2 = M_.facet_corners.adjacent_facet(c2);
        index_t v3 = M_.facet_corners.vertex(c3);

        M_.facet_corners.set_vertex(c1,v3);
        M_.facet_corners.set_adjacent_facet(c1,f2);
        M_.facet_corners.set_adjacent_facet(c2,f1);
        M_.facet_corners.set_vertex(c3,v1);

        // Re-insert t into the additional combinatorial data structure.
        insert(t);
    }

    /**
         * \brief Inserts a triangle of the mesh into the data structures
         *  used for topology checks.
         * \param[in] t index of the triangles to be inserted
         * \pre \p t is a valid triangle index in the mesh
     */
    void insert(index_t t) {
        for(
            index_t c=M_.facets.corners_begin(t);
            c<M_.facets.corners_end(t); ++c
        ) {
            index_t v = M_.facet_corners.vertex(c);
            if(v2c_[v] == NO_CORNER) {
                v2c_[v] = c;
                next_c_around_v_[c] = c;
            } else {
                next_c_around_v_[c] = next_c_around_v_[v2c_[v]];
                next_c_around_v_[v2c_[v]] = c;
            }
        }
    }

    /**
         * \brief Removes a triangle of the mesh from the data structures
         *  used for topology/combinatorial checks.
         * \param[in] t index of the triangles to be removed
         * \param[in] disconnect if true, connections from the neighbors
         *  to t are set to -1 (facet_corners.adjacent_facet).
         * \pre \p t is a valid triangle index in the mesh
     */
    void remove(index_t t, bool disconnect=true) {
        if(disconnect) {
            for(
                index_t c=M_.facets.corners_begin(t);
                c<M_.facets.corners_end(t); ++c
            ) {

                // Disconnect facet-facet link that point to t
                index_t t2 = M_.facet_corners.adjacent_facet(c);
                if(t2 != NO_FACET) {
                    for(
                        index_t c2=M_.facets.corners_begin(index_t(t2));
                        c2<M_.facets.corners_end(index_t(t2));
                        ++c2
                    ) {
                        if(
                            M_.facet_corners.adjacent_facet(c2) == t
                        ) {
                            M_.facet_corners.set_adjacent_facet(
                                c2,NO_FACET
                            );
                        }
                    }
                }
            }
        }


        for(
            index_t c=M_.facets.corners_begin(t);
            c<M_.facets.corners_end(t); ++c
        ) {
            // Remove t from combinatorial data structures
            index_t v = M_.facet_corners.vertex(c);
            if(next_c_around_v_[c] == c) {
                v2c_[v] = NO_CORNER;
            } else {
                index_t c_pred = next_c_around_v_[c];
                while(next_c_around_v_[c_pred] != c) {
                    c_pred = next_c_around_v_[c_pred];
                }
                next_c_around_v_[c_pred] = next_c_around_v_[c];
                v2c_[v] = c_pred;
            }
        }
    }

    /**
         * \brief Gets the number of triangles incident
         *  to a vertex.
         * \param[in] v index of the vertex
         * \return the number of triangles incident to \p v
     */
    index_t nb_incident_triangles(index_t v) const {
        index_t result = 0;
        index_t c = v2c_[v];
        do {
            ++result;
            c = next_c_around_v_[c];
        } while(c != v2c_[v]);
        return result;
    }

    /**
         * \brief Tests whether a given vertex is non-manifold
         *  by excess.
         * \details A vertex is non-manifold by-excess if its
         *  set of incident triangles contains a closed loop
         *  of triangles and additional triangles.
         * \param[in] v index of the vertex to be tested
         * \retval true if \p v is non-manifold by excess
         * \retval false otherwise
     */
    bool vertex_is_non_manifold_by_excess(index_t v, bool& moebius) {
        index_t nb_v_neighbors = nb_incident_triangles(v);
        index_t c = v2c_[v];
        do {
            index_t loop_size=0;
            index_t c_cur = c ;
            do {
                ++loop_size;
                if(c_cur == NO_CORNER) {
                    break;
                }
                if(loop_size > 100) {
                    // Probably Moebious strip or something...
                    moebius = true;
                    break;
                }
                c_cur = next_around_vertex_unoriented(v,c_cur);
            } while(c_cur != c);

            if(c_cur == c && loop_size < nb_v_neighbors) {
                return true;
            }
            c = next_c_around_v_[c];
        } while(c != v2c_[v]);

        return false;
    }

    /**
         * \brief Gets the next corner around a vertex from a given
         *  corner.
         * \details This function works even for a mesh that has triangles
         *  that are not coherently oriented. In other words, for two
         *  corners c1, c2, if we have:
         *   - v1 = facet_corners.vertex(c1)
         *   - v2 = facet_corners.vertex(
         *       c1,facets.next_corner_around_facet(c2f(c1),c1)
         *   )
         *   - w1 = facet_corners.vertex(c2)
         *   - w2 = facet_corners.vertex(
         *         c2,facets.next_corner_around_facet(c2f(c2),c2)
         *   )
         *  then we can have:
         *   - v1=w2 and v2=w1 (as usual) or:
         *   - v1=v2 and w1=w2 ('inverted' configuration)
         * \param[in] v the vertex
         * \param[in] c1 a corner incident to \p v or pointing to \p v
         * \return another corner incident to the \p v
     */
    index_t next_around_vertex_unoriented(
        index_t v, index_t c1
    ) const {
        index_t f1 = c2f(c1);
        index_t v1 = M_.facet_corners.vertex(c1);
        index_t v2 = M_.facet_corners.vertex(
            M_.facets.next_corner_around_facet(f1,c1)
        );

        geo_debug_assert(v1 == v || v2 == v);

        index_t f2 = M_.facet_corners.adjacent_facet(c1);
        if(f2 != NO_FACET) {
            for(
                index_t c2 = M_.facets.corners_begin(f2);
                c2 < M_.facets.corners_end(f2);
                ++c2
            ) {
                index_t w1 = M_.facet_corners.vertex(c2);
                index_t w2 = M_.facet_corners.vertex(
                    M_.facets.next_corner_around_facet(f2,c2)
                );
                if(
                    (v1 == w1 && v2 == w2) ||
                    (v1 == w2 && v2 == w1)
                ) {
                    if(w2 == v) {
                        return M_.facets.next_corner_around_facet(f2,c2);
                    } else {
                        geo_debug_assert(w1 == v);
                        return M_.facets.prev_corner_around_facet(f2,c2);
                    }
                }
            }
        }
        return NO_CORNER;
    }

    /**
         * \brief Gets the three corners adjacent to a triangle.
         * \details This function works even for a mesh that has triangles
         *  that are not coherently oriented. In other words, for two
         *  corners c1, c2, if we have:
         *   - v1 = facet_corners.vertex(c1)
         *   - v2 = facet_corners.vertex(
         *       c1,facets.next_corner_around_facet(c2f(c1),c1))
         *   - w1 = facet_corners.vertex(c2)
         *   - w2 = facet_corners.vertex(
         *       c2,facets.next_corner_around_facet(c2f(c2),c2))
         *  then c1 and c2 are adjacent if we have:
         *   - v1=w2 and v2=w1 (as usual) or:
         *   - v1=v2 and w1=w2 ('inverted' configuration)
         * \param[in] t1 index of the triangle
         * \param[out] adj_c index of the adjacent corners
         *  (array of 3 integers). Each entry contains a valid corner index
         *  or NO_CORNER if the corresponding edge is on the border.
         * \retval true if the three edges are manifold
         * \retval false otherwise (and then \p adj_c contains undefined
         *  values).
     */
    bool get_adjacent_corners(index_t t1, index_t* adj_c) {
        for(
            index_t c1 = M_.facets.corners_begin(t1);
            c1<M_.facets.corners_end(t1); ++c1
        ) {
            index_t v2 = M_.facet_corners.vertex(
                M_.facets.next_corner_around_facet(t1,c1)
            );

            *adj_c = NO_CORNER;

            // Traverse the circular incident edge list
            index_t c2=next_c_around_v_[c1];
            while(c2 != c1) {
                index_t t2 = c2f(c2);
                index_t c3 = M_.facets.prev_corner_around_facet(t2,c2);
                index_t v3 = M_.facet_corners.vertex(c3);
                if(v3 == v2) {
                    // Found an adjacent edge
                    if(*adj_c == NO_CORNER) {
                        *adj_c = c3;
                        geo_debug_assert(c3 != c1);
                    } else {
                        // If there was already an adjacent edge,
                        // then this is a non-manifold configuration
                        return false;
                    }
                }

                // Check with the other (wrong) orientation
                c3 = M_.facets.next_corner_around_facet(t2,c2);
                v3 = M_.facet_corners.vertex(c3);
                if(v3 == v2) {
                    // Found an adjacent edge
                    if(*adj_c == NO_CORNER) {
                        *adj_c = c2;
                        geo_debug_assert(c2 != c1);
                    } else {
                        // If there was already an adjacent edge,
                        // then this is a non-manifold configuration
                        return false;
                    }
                }
                c2 = next_c_around_v_[c2];
            }
            ++adj_c;
        }
        return true;
    }

    /**
         * \brief Tentatively connect a triangle of the mesh with its
         *   neighbors.
         * \details This function is independent of triangles orientations,
         *  see get_adjacent_corners().
         * \param[in] t index of the triangle to be connected
         * \param[in] adj_c an array of three integers that indicate
         *  for each corner of the triangle the index of the adjacent
         *  corner or NO_CORNER if the corner is on the border.
     */
    void connect_adjacent_corners(index_t t, index_t* adj_c) {
        for(index_t i=0; i<3; ++i) {
            if(adj_c[i] != NO_CORNER) {
                index_t c = M_.facets.corners_begin(t)+i;
                M_.facet_corners.set_adjacent_facet(c, c2f(adj_c[i]));
                M_.facet_corners.set_adjacent_facet(adj_c[i], t);
            }
        }
    }

    /**
         * \brief Tentatively connect a triangle of the mesh with its
         *   neighbors.
         * \details This function is independent of triangles orientations,
         *  see get_adjacent_corners().
         * \param[in] t index of the triangle to be connected
         * \retval false if the connection would have created non-manifold
         *  edges
         * \retval true otherwise
     */
    bool connect(index_t t) {
        index_t adj_c[3];
        if(!get_adjacent_corners(t,adj_c)) {
            return false;
        }
        connect_adjacent_corners(t, adj_c);
        return true;
    }

    /**
         * \brief Gets a facet index by corner index.
         * \details for a triangulated mesh, indexing is
         *  implicit, and we do not need to store a c2f array.
         * \param[in] c corner index
         * \return the index of the facet incident to c
     */
    index_t c2f(index_t c) const {
        geo_debug_assert(c != NO_CORNER);
        geo_debug_assert(c < M_.facet_corners.nb());
        return c/3;
    }


    /**
         * \brief Tests whether two triangles have the
         *  same orientation.
         * \param[in] t1 first triangle
         * \param[in] t2 second triangle
         * \retval true if \p t1 and \p t2 have the same
         *  orientation
         * \retval false otherwise
         * \pre \p t1 and \p t2 share an edge
     */
    bool triangles_have_same_orientation(
        index_t t1,
        index_t t2
    ) {
        index_t c1 = M_.facets.corners_begin(t1);
        index_t i1 = M_.facet_corners.vertex(c1);
        index_t j1 = M_.facet_corners.vertex(c1+1);
        index_t k1 = M_.facet_corners.vertex(c1+2);

        index_t c2 = M_.facets.corners_begin(t2);
        index_t i2 = M_.facet_corners.vertex(c2);
        index_t j2 = M_.facet_corners.vertex(c2+1);
        index_t k2 = M_.facet_corners.vertex(c2+2);

        if(
            (i1==i2 && j1==j2) ||
            (i1==k2 && j1==i2) ||
            (i1==j2 && j1==k2) ||
            (k1==k2 && i1==i2) ||
            (k1==j2 && i1==k2) ||
            (k1==i2 && i1==j2) ||
            (j1==j2 && k1==k2) ||
            (j1==i2 && k1==j2) ||
            (j1==k2 && k1==i2)
        ) {
            return false;
        }

        return true;
    }


    /**
         * \brief Tests whether the normals of two triangles that
         *  share an edge 'agree', i.e. whether they do not form
         *  a too sharp angle.
         * \param[in] t1 index of the first triangle
         * \param[in] t2 index of the second triangle
         * \retval true if the normals of both triangles do not
         *  point in opposite directions
         * \retval false otherwise
         * \pre the two triangles are incident to the same edge
         *  (they have two vertices in common)
     */
    bool triangles_normals_agree(
        index_t t1,
        index_t t2
    ) const {
        const vec3* points =
            reinterpret_cast<const vec3*>(M_.vertices.point_ptr(0));

        index_t c1 = M_.facets.corners_begin(t1);
        index_t i1 = M_.facet_corners.vertex(c1);
        index_t j1 = M_.facet_corners.vertex(c1+1);
        index_t k1 = M_.facet_corners.vertex(c1+2);

        index_t c2 = M_.facets.corners_begin(t2);
        index_t i2 = M_.facet_corners.vertex(c2);
        index_t j2 = M_.facet_corners.vertex(c2+1);
        index_t k2 = M_.facet_corners.vertex(c2+2);

        vec3 n1 = normalize(
            cross(
                points[j1] - points[i1],
                points[k1] - points[i1]
                )
        );

        vec3 n2 = normalize(
            cross(
                points[j2] - points[i2],
                points[k2] - points[i2]
                )
        );

        double d = dot(n1,n2);
        // Test for combinatorial orientation,
        // if t1 and t2 have opposite orientation,
        // then we flip one of the normals (i.e.,
        // we simply change the sign of the dot product).
        if(
            (i1==i2 && j1==j2) ||
            (i1==k2 && j1==i2) ||
            (i1==j2 && j1==k2) ||
            (k1==k2 && i1==i2) ||
            (k1==j2 && i1==k2) ||
            (k1==i2 && i1==j2) ||
            (j1==j2 && k1==k2) ||
            (j1==i2 && k1==j2) ||
            (j1==k2 && k1==i2)
        ) {
            d = -d;
        }
        return (d > -0.8);
    }

    /**
         * \brief Merges two connected components.
         * \details The connected component incident to \p t
         *  is replaced with \p comp2.
         * \param [in] t index of a triangle incident
         *  to the first connected component
         * \param [in] comp2 index of the second connected
         *  component
         * \param [in] flip if true, flip the triangles
         * \pre At least one of the triangles adjacent to
         *  \p t (directly or not) is incident to
         *  component \p comp2
     */
    void merge_connected_component(
        index_t t,
        index_t comp2,
        bool flip
    ) {
        geo_assert(comp2 != cnx_[t]);

        std::stack<index_t> S;
        index_t comp1 = cnx_[t];


        cnx_[t] = comp2;
        --cnx_size_[comp1];
        ++cnx_size_[comp2];
        if(flip) {
            flip_triangle(t);
        }
        S.push(t);
        while(!S.empty()) {
            index_t t1 = S.top();
            S.pop();
            for(
                index_t c = M_.facets.corners_begin(t1);
                c < M_.facets.corners_end(t1); ++c
            ) {
                index_t t2 = M_.facet_corners.adjacent_facet(c);
                if(t2 != NO_FACET && cnx_[t2] == comp1) {
                    cnx_[t2] = comp2;
                    --cnx_size_[comp1];
                    ++cnx_size_[comp2];
                    if(flip) {
                        flip_triangle(t2);
                    }
                    S.push(t2);
                }
            }
        }
        geo_assert(cnx_size_[comp1] == 0);
    }

    /**
         * \brief Initializes the date structures
         *  that represent the connected components.
         * \details This function computes cnx_ and
         *  cnx_size_. The array cnx_[f] gives for each
         *  facet f the index of the connected component
         *  that contains f, and the array cnx_size_[comp]
         *  gives for each connected component comp the
         *  number of facets in comp.
     */
    void init_connected_components() {
        cnx_.assign(M_.facets.nb(), index_t(NO_CNX));
        cnx_size_.clear();
        for(index_t t=0; t<M_.facets.nb(); ++t) {
            if(cnx_[t] == NO_CNX) {
                index_t cnx_id  = cnx_size_.size();
                index_t nb = 0;
                std::stack<index_t> S;
                S.push(t);
                cnx_[t] = cnx_id;
                ++nb;
                while(!S.empty()) {
                    index_t t2 = S.top();
                    S.pop();
                    for(
                        index_t c=M_.facets.corners_begin(t2);
                        c<M_.facets.corners_end(t2); ++c
                    ) {
                        index_t t3 = M_.facet_corners.adjacent_facet(c);
                        if(t3 != NO_FACET && cnx_[t3] != cnx_id) {
                            geo_assert(cnx_[t3] == NO_CNX);
                            cnx_[t3] = cnx_id;
                            ++nb;
                            S.push(t3);
                        }
                    }
                }
                cnx_size_.push_back(nb);
            }
        }
        Logger::out("Co3Ne")
            << "Found " << cnx_size_.size() << " connected components"
            << std::endl;
    }

  private:
    Mesh& M_;

    /**
         * \brief For each corner, next_c_around_v_[c]
         * chains the circular list of corners
         * incident to the same corner as c.
     */
    vector<index_t> next_c_around_v_;

    /**
         * \brief For each vertex v, v2c_[v] contains a
         *  corner incident to v, or NO_VERTEX if v is
         *  isolated.
     */
    vector<index_t> v2c_;


    /**
         * \brief For each triangle t, cnx_[t] contains
         *  the index of the connected component of the
         *  mesh incident to t.
     */
    vector<index_t> cnx_;

    /**
         * \brief For each connected component C,
         *  cnx_size_[C] contains the number of
         *  facets in C.
     */
    vector<index_t> cnx_size_;

    /**
         * \brief In strict mode, each inserted triangle
         *  is checked for non-manifold configuration.
         *  In non-strict mode, only T2 and T1 triangles are
         *  tested (those seen from only 2 or only 1 Voronoi
         *  cell), T3 triangles are inserted without test.
     */
    bool strict_;
};
}
#endif // LPCVT_RECON_CO3NEMANIFOLDEXTRACTION_H
