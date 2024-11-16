//
// Created by hongbo on 15/11/24.
//

#ifndef LPCVT_RECON_CO3NERVD_H
#define LPCVT_RECON_CO3NERVD_H
#include <geogram/basic/geometry.h>
#include <geogram/points/nn_search.h>
#include <geogram/mesh/mesh.h>
#include "../CVT_UTILS/Utils.h"
namespace NMCVT_PC{
    using namespace GEO;

class Co3NeRestrictedVoronoiDiagram {
  public:
    /**
     * \brief number of elements in sine/cosine table.
     */
    static const index_t sincos_nb = 10;


    /**
         * \brief Stores a 3D point and the combinatorial information
         *  (index of the adjacent seed). The combinatorial information
         *  is used to reconstruct the triangles at the end of the
         *  algorithm.
     */
    class Vertex {
      public:
        /**
             * \brief Constructs a new uninitialized Vertex.
         */
        Vertex() {
        }

        /**
             * \brief Constructs a Vertex from a 3d point.
         */
        Vertex(const vec3& v) :
                                point_(v),
                                adjacent_seed_(-1) {
        }

        /**
             * \brief Gets the 3d point associated with this vertex.
             * \return a const reference to the 3d point
         */
        const vec3& point() const {
            return point_;
        }

        /**
             * \brief Gets the 3d point associated with this vertex.
             * \return a const reference to the 3d point
         */
        vec3& point() {
            return point_;
        }

        /**
             * \brief Gets the index of the adjacent seed associated with
             *  this vertex.
             * \details Each vertex stores combinatorial information, i.e.
             *  the index of the adjacent Voronoi seed accros the edge
             *  starting from this vertex
             * \return the index of the adjacent Voronoi seed
         */
        signed_index_t adjacent_seed() const {
            return adjacent_seed_;
        }

        /**
             * \brief Sets the index of the adjacent seed associated with
             *  this vertex.
             * \details Each vertex stores combinatorial information, i.e.
             *  the index of the adjacent Voronoi seed accros the edge
             *  starting from this vertex
             * \param[in] x the index of the adjacent Voronoi seed
         */
        void set_adjacent_seed(signed_index_t x) {
            adjacent_seed_ = x;
        }

      private:
        vec3 point_;
        signed_index_t adjacent_seed_;
    };

    /**
         * \brief Internal representation of the polygons, that represent
         *  the intersection between the disks and the Voronoi cells.
     */
    class Polygon {
      public:
        /**
             * \brief Creates a new uninitialized polygon with a given
             *  number of vertices.
             * \param[in] size number of vertices
         */
        Polygon(index_t size) :
                                vertices_(size) {
        }

        /**
             * \brief Gets the number of vertices.
             * \return the number of vertices of this Polygon
         */
        index_t nb_vertices() const {
            return vertices_.size();
        }

        /**
             * \brief Adds a new vertex to this Polygon.
             * \param[in] v the vertex to be added.
         */
        void add_vertex(const Vertex& v) {
            vertices_.push_back(v);
        }

        /**
             * \brief Gets a Vertex by its index.
             * \param[in] i the index of the Vertex
             * \return a reference to the Vertex
         */
        Vertex& vertex(index_t i) {
            return vertices_[i];
        }

        /**
             * \brief Gets a Vertex by its index.
             * \param[in] i the index of the Vertex
             * \return a const reference to the Vertex
         */
        const Vertex& vertex(index_t i) const {
            return vertices_[i];
        }

        /**
             * \brief Gets the index of the next vertex around
             *  the polygon.
             * \param[in] i index of the vertex
             * \return index of the next vertex (successor of \p i)
             *  around the Polygon.
         */
        index_t next_vertex(index_t i) const {
            return (i == nb_vertices() - 1) ? 0 : i + 1;
        }

        /**
             * \brief Removes all the vertices.
         */
        void clear() {
            vertices_.resize(0);
        }

        /**
             * \brief Swaps the vertices of this Polygon with
             *  the vertices of another polygon.
             * \param[in] P the other polygon
         */
        void swap(Polygon& P) {
            vertices_.swap(P.vertices_);
        }

      private:
        vector<Vertex> vertices_;
    };

    /**
         * \brief Constructs a new uninitialized Co3NeRestrictedVoronoiDiagram.
     */
    Co3NeRestrictedVoronoiDiagram() :
                                      nb_points_(0),
                                      p_(nullptr),
                                      p_stride_(0),
                                      n_(nullptr),
                                      n_stride_(0),
                                      radius_(0.0),
                                      NN_(NearestNeighborSearch::create(3)),
                                      sqROS_(0.0),
                                      nb_neighbors_(0)
    {
    }

    /**
         * \brief Co3NeRestrictedVoronoiDiagram destructor.
     */
    ~Co3NeRestrictedVoronoiDiagram() {
        clear();
    }


    /**
         * \brief Sets or resets exact mode for nearest neighbor search
         * (default is exact).
         * \details Nearest neighbor search can be exact or approximate.
         * Note that approximate mode cannot be used for the
         * final reconstruction phase (that needs exact combinatorics),
         * but it may speedup the smoothing phase.
         * \param[in] x if set, nearest neighbors search are exact, else they
         *  are approximate
     */
    void set_exact(bool x) {
        NN_->set_exact(x);
    }

    /**
         * \brief Clears this Co3NeRestrictedVoronoiDiagram.
     */
    void clear() {
        NN_.reset();
        nb_points_ = 0;
        p_ = nullptr;
        p_stride_ = 0;
        n_ = nullptr;
        n_stride_ = 0;
        nb_neighbors_ = 0;
    }

    /**
         * \brief Initializes this Co3NeRestrictedVoronoiDiagram from a
         *  pointset stored in a mesh.
         * \details If the mesh \p M has normals, then they are used.
         * \param[in] M the pointset
     */
    void init(Mesh& M) {
        geo_assert(M.vertices.dimension() >= 3);
        geo_assert(M.vertices.dimension() == 3 || NN_->stride_supported());
        double* normals_pointer = nullptr;
        {
            Attribute<double> normal;
            normal.bind_if_is_defined(M.vertices.attributes(), "normal");
            if(normal.is_bound() && normal.dimension() == 3) {
                normals_pointer = &normal[0];
            }
        }

        if(normals_pointer == nullptr) {
            init(
                M.vertices.nb(),
                M.vertices.point_ptr(0), M.vertices.dimension(),
                nullptr, 0
            );
        } else {
            init(
                M.vertices.nb(),
                M.vertices.point_ptr(0), M.vertices.dimension(),
                normals_pointer, 3
            );
        }
    }

    /**
         * \brief Initializes this Co3NeRestrictedVoronoiDiagram from an
         *  array of points and an array of normals.
         * \param[in] nb_points_in number of points
         * \param[in] p pointer to the coordinates of the points
         * \param[in] p_stride number of doubles between two consecutive points
         * \param[in] n pointer to the coordinates of the normals
         * \param[in] n_stride number of doubles between two consecutive normals
     */
    void init(
        index_t nb_points_in,
        double* p, index_t p_stride,
        double* n, index_t n_stride
    ) {
        nb_points_ = nb_points_in;
        p_ = p;
        p_stride_ = p_stride;
        n_ = n;
        n_stride_ = n_stride;
        NN_->set_points(nb_points(), p_, p_stride_);
    }

    /**
         * \brief Reconstructs the nearest neighbors search data
         * structure.
         * \details This function needs to be called whenever
         * the point set changes.
     */
    void update() {
        init(nb_points_, p_, p_stride_, n_, n_stride_);
    }

    /**
         * \brief Sets the radius of the circles used to determine
         *  points adjacencies.
         * \param[in] r the radius of the circles
     */
    void set_circles_radius(double r) {
        radius_ = r;
        sqROS_ = 4.0 * radius_ * radius_;  // squared radius of security
        // when a neighbor is further away than ROS, then it cannot
        // clip a circle of radius r
    }

    /**
         * \brief Gets the number of points.
         * \return the number of points
     */
    index_t nb_points() const {
        return nb_points_;
    }

    /**
         * \brief Gets a point by its index.
         * \param[in] i index of the point
         * \return a const reference to the point
     */
    const vec3& point(index_t i) const {
        geo_debug_assert(i < nb_points());
        return *(vec3*) (p_ + i * p_stride_);
        // Yes I know, this is a bit ugly...
    }

    /**
         * \brief Gets a normal by point index.
         * \param[in] i index of the point
         * \return a const reference to the normal
         *  associated with the point
     */
    const vec3& normal(index_t i) const {
        geo_debug_assert(n_ != nullptr);
        geo_debug_assert(i < nb_points());
        return *(vec3*) (n_ + i * n_stride_);
        // Yes I know, this is a bit ugly...
    }

    /**
         * \brief Sets the normal associated wigth a point.
         * \param[in] i the index of the point
         * \param[in] N the normal
     */
    void set_normal(index_t i, const vec3& N) const {
        geo_debug_assert(n_ != nullptr);
        geo_debug_assert(i < nb_points());
        double* n = n_ + i * n_stride_;
        n[0] = N.x;
        n[1] = N.y;
        n[2] = N.z;
    }

    /**
         * \brief Computes the intersection between a polygon
         *  and the halfspace determined by the bisector
         *  of two points.
         * \param[in,out] Ping polygon to be clipped
         * \param[in,out] Pong a temporary work variable provided
         *  by the caller
         * \param[in] pi first extremity of the bisector
         * \param[in] pj second extremity of the bisector
         * \param[in] j index of the second extremity of the
         *  bisector (used to store the combinatorial information).
     */
    static void clip_polygon_by_bisector(
        Polygon& Ping, Polygon& Pong,
        const vec3& pi, const vec3& pj, index_t j
    ) {
        if(Ping.nb_vertices() == 0) {
            return;
        }
        Pong.clear();

        vec3 n(
            pi.x - pj.x,
            pi.y - pj.y,
            pi.z - pj.z
        );

        // Compute d = n . m, where n is the
        // normal vector of the bisector [pi,pj]
        // and m twice the middle point of the bisector.
        double d =
            n.x * (pi.x + pj.x) +
            n.y * (pi.y + pj.y) +
            n.z * (pi.z + pj.z);

        // The predecessor of the first vertex is the last vertex
        index_t prev_k = Ping.nb_vertices() - 1;
        const Vertex* prev_vk = &(Ping.vertex(prev_k));

        // We compute:
        //    prev_l = prev_vk . n
        double prev_l = dot(prev_vk->point(), n);

        // We compute:
        //    side1(pi,pj,q) = sign(2*q.n - n.m) = sign(2*l - d)
        Sign prev_status = geo_sgn(2.0 * prev_l - d);

        for(index_t k = 0; k < Ping.nb_vertices(); k++) {
            const Vertex* vk = &(Ping.vertex(k));

            // We compute: l = vk . n
            double l = dot(vk->point(), n);

            // We compute:
            //   side1(pi,pj,q) = sign(2*q.n - n.m) = sign(2*l - d)
            Sign status = geo_sgn(2.0 * l - d);

            // If status of edge extremities differ,
            // then there is an intersection.
            if(status != prev_status && (prev_status != 0)) {

                // Compute lambda1 and lambda2, the
                // barycentric coordinates of the intersection I
                // in the segment [prev_vk vk]
                // Note that d and l (used for the predicates)
                // are reused here.
                double denom = 2.0 * (prev_l - l);
                double lambda1, lambda2;

                // Shit happens ! [Forrest Gump]
                if(::fabs(denom) < 1e-20) {
                    lambda1 = 0.5;
                    lambda2 = 0.5;
                } else {
                    lambda1 = (d - 2.0 * l) / denom;
                    // Note: lambda2 is also given
                    // by (2.0*l2-d)/denom
                    // (but 1.0 - lambda1 is a bit
                    //  faster to compute...)
                    lambda2 = 1.0 - lambda1;
                }
                Vertex V;
                V.point().x =
                    lambda1 * prev_vk->point().x + lambda2 * vk->point().x;
                V.point().y =
                    lambda1 * prev_vk->point().y + lambda2 * vk->point().y;
                V.point().z =
                    lambda1 * prev_vk->point().z + lambda2 * vk->point().z;
                if(status > 0) {
                    V.set_adjacent_seed(prev_vk->adjacent_seed());
                } else {
                    V.set_adjacent_seed(signed_index_t(j));
                }
                Pong.add_vertex(V);
            }
            if(status > 0) {
                Pong.add_vertex(*vk);
            }
            prev_vk = vk;
            prev_status = status;
            prev_k = k;
            prev_l = l;
        }
        Ping.swap(Pong);
    }

    /**
         * \brief Computes the squared maximum distance between a point
         *  and the vertices of a polygon.
         * \param[in] p the point
         * \param[in] P the polygon
         * \return the maximum squared distance between \p p and the vertices
         *  of \p P
     */
    static double squared_radius(const vec3& p, const Polygon& P) {
        double result = 0.0;
        for(index_t i = 0; i < P.nb_vertices(); i++) {
            result = std::max(result, distance2(p, P.vertex(i).point()));
        }
        return result;
    }

    /**
         * \brief Computes a polygon that approximates a disk centered
         *  at a point and orthogonal to its normal vector.
         * \param[in] i index of the point
         * \param[out] P an approximation of the circle centered
         *  at point \p i with normal vector \p N. The radius is
         *  defined by set_circles_radius().
         * \param[in] N normal vector
     */
    void get_circle(index_t i, Polygon& P, const vec3& N) const {
        P.clear();
        const vec3& pi = point(i);
        vec3 U = Geom::perpendicular(N);
        U = normalize(U);
        vec3 V = cross(N, U);
        V = normalize(V);
        // We use a table for sine and cosine for speeding up things
        // a little bit (especially on some cell phones / handheld devices
        // that do not have a FPU).
        /*
                    const index_t nb = 10;
                    for(index_t k=0; k<nb; ++k) {
                        double alpha = 2.0 * M_PI * double(k) / double(nb - 1);
                        double s = sin(alpha);
                        double c = cos(alpha);
                        vec3 p = pi + c * radius_ * U + s * radius_ * V;
                        P.add_vertex(p);
                    }
        */

        for(index_t k = 0; k < sincos_nb; ++k) {
            double s = sincos_table[k][0];
            double c = sincos_table[k][1];
            vec3 p = pi + c * radius_ * U + s * radius_ * V;
            P.add_vertex(p);
        }

    }

    /**
         * \brief Nearest neighbor search
         * \param[in] i index of the query point
         * \param[out] neigh array of nb signed_index_t
         * \param[out] sq_dist array of nb doubles
         * \param[in] nb number of neighbors to be searched
     */
    void get_neighbors(
        index_t i,
        index_t* neigh,
        double* sq_dist,
        index_t nb
    ) const {
        return NN_->get_nearest_neighbors(
            nb, i, neigh, sq_dist
        );
    }

    /**
         * \brief Nearest neighbor search
         * \param[in] i index of the query point
         * \param[out] neigh vector of signed_index_t
         * \param[out] sq_dist array of nb doubles
         * \param[in] nb number of neighbors to be searched
     */
    void get_neighbors(
        index_t i,
        vector<index_t>& neigh,
        vector<double>& sq_dist,
        index_t nb
    ) const {
        neigh.resize(nb);
        sq_dist.resize(nb);
        get_neighbors(i, neigh.data(), sq_dist.data(), nb);
    }

    /**
         * \brief Computes a Restricted Voronoi Cell (RVC), i.e.
         *  the intersection between a disk and the Voronoi cell
         *  of a point.
         * \details The temporary work variables provided by the caller
         *  make it possible to reuse memory accros multiple calls to this
         *  function and thus avoid multiple dynamic memory allocations.
         * \param[in] i index of the point that determines the Voronoi cell.
         * \param[out] P result
         * \param[in] Q work temporary variable provided by caller
         * \param[in] neighbor work temporary variable provided by caller
         * \param[in] squared_dist work temporary variable provided by caller
     */
    void get_RVC(
        index_t i, Polygon& P,
        Polygon& Q,
        vector<index_t>& neighbor,
        vector<double>& squared_dist
    ) const {
        neighbor.resize(0);
        squared_dist.resize(0);
        get_RVC(i, normal(i), P, Q, neighbor, squared_dist);
    }

    /**
         * \brief Computes a Restricted Voronoi Cell (RVC), i.e.
         *  the intersection between a disk and the Voronoi cell
         *  of a point.
         * \details The temporary work variables provided by the caller
         *  make it possible to reuse memory accros multiple calls to this
         *  function and thus avoid multiple dynamic memory allocations.
         * \param[in] i index of the point that determines the Voronoi cell.
         * \param[in] N normal vector at point \p i
         * \param[out] P result
         * \param[in] Q work temporary variable, provided by caller
         * \param[in] neighbor initial neighbor indices
         *  if size is not zero, contains (previously computed)
         *  neighbor indices.
         * \param[in] squared_dist initial neighbor squared distances
         *  if size is not zero, contains (previously computed)
         *  neighbor squared distances.
     */
    void get_RVC(
        index_t i, const vec3& N, Polygon& P,
        Polygon& Q,
        vector<index_t>& neighbor,
        vector<double>& squared_dist
    ) const {
        get_circle(i, P, N);

        index_t nb_neigh = std::min(index_t(nb_points() - 1), index_t(20));
        index_t jj = 0;

        // just in case, limit to 1000 neighbors.
        index_t max_neigh = std::min(index_t(1000), nb_points() - 1);

        while(nb_neigh < max_neigh) {
            if(P.nb_vertices() < 3) {
                return;
            }
            if(neighbor.size() < nb_neigh) {
                get_neighbors(i, neighbor, squared_dist, nb_neigh);
            }
            while(jj < nb_neigh && squared_dist[jj] < 1e-30) {
                jj++;
            }
            while(jj < nb_neigh) {
                if(squared_dist[jj] > sqROS_) {
                    return;
                }
                index_t j = neighbor[jj];
                double Rk = squared_radius(point(i), P);
                if(squared_dist[jj] > 4.0 * Rk) {
                    return;
                }
                clip_polygon_by_bisector(P, Q, point(i), point(j), j);
                jj++;
            }
            if(nb_neigh > 3) {
                nb_neigh += nb_neigh / 3;
            } else {
                nb_neigh++;
            }
            nb_neigh = std::min(nb_neigh, nb_points()-1);
        }
    }

    /**
         * \brief Gets the number of neighbors, used for nearest neighbors
         *  queries.
         * \return the number of neighbors
     */
    index_t nb_neighbors() const {
        return std::min(nb_neighbors_,nb_points()-1);
    }

    /**
         * \brief Sets the number of neighbors, used for nearest neighbors
         *  queries.
         * \param[in] x the number of neighbors
     */
    void set_nb_neighbors(index_t x) {
        nb_neighbors_ = x;
    }

  private:
    friend class Co3Ne;

    index_t nb_points_;
    double* p_;
    index_t p_stride_;
    double* n_;
    index_t n_stride_;
    double radius_;

    NearestNeighborSearch_var NN_;

    double sqROS_;
    index_t nb_neighbors_;
};

}

#endif // LPCVT_RECON_CO3NERVD_H
