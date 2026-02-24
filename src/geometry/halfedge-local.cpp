
#include "halfedge.h"

#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <iostream>

/******************************************************************
*********************** Local Operations **************************
******************************************************************/

/* Note on local operation return types:

    The local operations all return a std::optional<T> type. This is used so that your
    implementation can signify that it cannot perform an operation (i.e., because
    the resulting mesh does not have a valid representation).

    An optional can have two values: std::nullopt, or a value of the type it is
    parameterized on. In this way, it's similar to a pointer, but has two advantages:
    the value it holds need not be allocated elsewhere, and it provides an API that
    forces the user to check if it is null before using the value.

    In your implementation, if you have successfully performed the operation, you can
    simply return the required reference:

            ... collapse the edge ...
            return collapsed_vertex_ref;

    And if you wish to deny the operation, you can return the null optional:

            return std::nullopt;

    Note that the stubs below all reject their duties by returning the null optional.
*/


/*
 * add_face: add a standalone face to the mesh
 *  sides: number of sides
 *  radius: distance from vertices to origin
 *
 * We provide this method as an example of how to make new halfedge mesh geometry.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::add_face(uint32_t sides, float radius) {
	//faces with fewer than three sides are invalid, so abort the operation:
	if (sides < 3) return std::nullopt;


	std::vector< VertexRef > face_vertices;
	//In order to make the first edge point in the +x direction, first vertex should
	// be at -90.0f - 0.5f * 360.0f / float(sides) degrees, so:
	float const start_angle = (-0.25f - 0.5f / float(sides)) * 2.0f * PI_F;
	for (uint32_t s = 0; s < sides; ++s) {
		float angle = float(s) / float(sides) * 2.0f * PI_F + start_angle;
		VertexRef v = emplace_vertex();
		v->position = radius * Vec3(std::cos(angle), std::sin(angle), 0.0f);
		face_vertices.emplace_back(v);
	}

	assert(face_vertices.size() == sides);

	//assemble the rest of the mesh parts:
	FaceRef face = emplace_face(false); //the face to return
	FaceRef boundary = emplace_face(true); //the boundary loop around the face

	std::vector< HalfedgeRef > face_halfedges; //will use later to set ->next pointers

	for (uint32_t s = 0; s < sides; ++s) {
		//will create elements for edge from a->b:
		VertexRef a = face_vertices[s];
		VertexRef b = face_vertices[(s+1)%sides];

		//h is the edge on face:
		HalfedgeRef h = emplace_halfedge();
		//t is the twin, lies on boundary:
		HalfedgeRef t = emplace_halfedge();
		//e is the edge corresponding to h,t:
		EdgeRef e = emplace_edge(false); //false: non-sharp

		//set element data to something reasonable:
		//(most ops will do this with interpolate_data(), but no data to interpolate here)
		h->corner_uv = a->position.xy() / (2.0f * radius) + 0.5f;
		h->corner_normal = Vec3(0.0f, 0.0f, 1.0f);
		t->corner_uv = b->position.xy() / (2.0f * radius) + 0.5f;
		t->corner_normal = Vec3(0.0f, 0.0f,-1.0f);

		//thing -> halfedge pointers:
		e->halfedge = h;
		a->halfedge = h;
		if (s == 0) face->halfedge = h;
		if (s + 1 == sides) boundary->halfedge = t;

		//halfedge -> thing pointers (except 'next' -- will set that later)
		h->twin = t;
		h->vertex = a;
		h->edge = e;
		h->face = face;

		t->twin = h;
		t->vertex = b;
		t->edge = e;
		t->face = boundary;

		face_halfedges.emplace_back(h);
	}

	assert(face_halfedges.size() == sides);

	for (uint32_t s = 0; s < sides; ++s) {
		face_halfedges[s]->next = face_halfedges[(s+1)%sides];
		face_halfedges[(s+1)%sides]->twin->next = face_halfedges[s]->twin;
	}

	return face;
}


/*
 * bisect_edge: split an edge without splitting the adjacent faces
 *  e: edge to split
 *
 * returns: added vertex
 *
 * We provide this as an example for how to implement local operations.
 * (and as a useful subroutine!)
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::bisect_edge(EdgeRef e) {
	// Phase 0: draw a picture
	//
	// before:
	//    ----h--->
	// v1 ----e--- v2
	//   <----t---
	//
	// after:
	//    --h->    --h2->
	// v1 --e-- vm --e2-- v2
	//    <-t2-    <--t--
	//

	// Phase 1: collect existing elements
	HalfedgeRef h = e->halfedge;
	HalfedgeRef t = h->twin;
	VertexRef v1 = h->vertex;
	VertexRef v2 = t->vertex;

	// Phase 2: Allocate new elements, set data
	VertexRef vm = emplace_vertex();
	vm->position = (v1->position + v2->position) / 2.0f;
	interpolate_data({v1, v2}, vm); //set bone_weights

	EdgeRef e2 = emplace_edge();
	e2->sharp = e->sharp; //copy sharpness flag

	HalfedgeRef h2 = emplace_halfedge();
	interpolate_data({h, h->next}, h2); //set corner_uv, corner_normal

	HalfedgeRef t2 = emplace_halfedge();
	interpolate_data({t, t->next}, t2); //set corner_uv, corner_normal

	// The following elements aren't necessary for the bisect_edge, but they are here to demonstrate phase 4
    FaceRef f_not_used = emplace_face();
    HalfedgeRef h_not_used = emplace_halfedge();

	// Phase 3: Reassign connectivity (careful about ordering so you don't overwrite values you may need later!)

	vm->halfedge = h2;

	e2->halfedge = h2;

	assert(e->halfedge == h); //unchanged

	//n.b. h remains on the same face so even if h->face->halfedge == h, no fixup needed (t, similarly)

	h2->twin = t;
	h2->next = h->next;
	h2->vertex = vm;
	h2->edge = e2;
	h2->face = h->face;

	t2->twin = h;
	t2->next = t->next;
	t2->vertex = vm;
	t2->edge = e;
	t2->face = t->face;
	
	h->twin = t2;
	h->next = h2;
	assert(h->vertex == v1); // unchanged
	assert(h->edge == e); // unchanged
	//h->face unchanged

	t->twin = h2;
	t->next = t2;
	assert(t->vertex == v2); // unchanged
	t->edge = e2;
	//t->face unchanged


	// Phase 4: Delete unused elements
    erase_face(f_not_used);
    erase_halfedge(h_not_used);

	// Phase 5: Return the correct iterator
	return vm;
}


/*
 * split_edge: split an edge and adjacent (non-boundary) faces
 *  e: edge to split
 *
 * returns: added vertex. vertex->halfedge should lie along e
 *
 * Note that when splitting the adjacent faces, the new edge
 * should connect to the vertex ccw from the ccw-most end of e
 * within the face.
 *
 * Do not split adjacent boundary faces.
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::split_edge(EdgeRef e) {
    HalfedgeRef h = e->halfedge;
    HalfedgeRef t = h->twin;

    VertexRef A = h->vertex;
    VertexRef B = t->vertex;
    FaceRef F_h = h->face;
    FaceRef F_t = t->face;

    HalfedgeRef hn = h->next;
    HalfedgeRef tn = t->next;

    HalfedgeRef h_prev = h;
    while (h_prev->next != h) {
        h_prev = h_prev->next;
    }
    HalfedgeRef t_prev = t;
    while (t_prev->next != t) {
        t_prev = t_prev->next;
    }

    VertexRef vm = emplace_vertex();
    vm->position = (A->position + B->position) / 2.0f;
    interpolate_data({A, B}, vm);

    EdgeRef e_new = emplace_edge();
    e_new->sharp = e->sharp;

    HalfedgeRef h2 = emplace_halfedge();
    HalfedgeRef t2 = emplace_halfedge();

    interpolate_data({h, hn}, h2);
    interpolate_data({t, tn}, t2);

    h->twin = t2;  
    t2->twin = h;
    h2->twin = t;   
    t->twin = h2;

    h->edge = e;       
    t2->edge = e;
    h2->edge= e_new;   
    t->edge = e_new;
    e->halfedge = h;
    e_new->halfedge = h2;

    h2->vertex = vm;
    t2->vertex = vm;

    vm->halfedge = t2;

    if (!F_h->boundary) {
        HalfedgeRef hnn = hn->next;
        
        FaceRef F_h_new = emplace_face(false);
        EdgeRef e_diag_h = emplace_edge(false);
        HalfedgeRef hd = emplace_halfedge();
        HalfedgeRef hd_t = emplace_halfedge();

        interpolate_data({h, hn}, hd);
        interpolate_data({hn, hnn}, hd_t);

        hd->twin = hd_t;  
        hd_t->twin = hd;
        hd->edge = e_diag_h;  
        hd_t->edge = e_diag_h;
        e_diag_h->halfedge = hd;

        hd->vertex = vm;          
        hd->face = F_h;
        hd_t->vertex = hnn->vertex;  
        hd_t->face = F_h_new;

        h->next = hd;
        hd->next = hnn;

        h2->face = F_h_new;
        h2->next = hn;
        hn->next = hd_t;
        hd_t->next = h2;
        hn->face = F_h_new;

        F_h->halfedge = h;
        F_h_new->halfedge = h2;
    } 
    else {
        h2->face = F_h;
        h->next = h2;
        h2->next = hn; 
        F_h->halfedge = h;
    }

    if (!F_t->boundary) {
        HalfedgeRef tnn = tn->next;

        FaceRef F_t_new = emplace_face(false);
        EdgeRef e_diag_t = emplace_edge(false);
        HalfedgeRef td = emplace_halfedge();
        HalfedgeRef td_t = emplace_halfedge();

        interpolate_data({t, tn}, td_t);
        interpolate_data({tn, tnn}, td);

        td->twin = td_t;  
        td_t->twin = td;
        td->edge = e_diag_t;  
        td_t->edge = e_diag_t;
        e_diag_t->halfedge = td;

        td->vertex = tnn->vertex;  
        td->face = F_t_new;
        td_t->vertex = vm;           
        td_t->face = F_t;

        t->next = td_t;
        td_t->next = tnn;

        t2->face = F_t_new;
        t2->next = tn;
        tn->next = td;
        td->next = t2;
        tn->face = F_t_new;

        F_t->halfedge = t;
        F_t_new->halfedge = t2;
    } 
    else {
        t2->face = F_t;
        t->next = t2;
        t2->next = tn;
        F_t->halfedge = t;
    }

    return vm;
}



/*
 * inset_vertex: divide a face into triangles by placing a vertex at f->center()
 *  f: the face to add the vertex to
 *
 * returns:
 *  std::nullopt if insetting a vertex would make mesh invalid
 *  the inset vertex otherwise
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::inset_vertex(FaceRef f) {
	// A2Lx4 (OPTIONAL): inset vertex
	
	(void)f;
    return std::nullopt;
}


/* [BEVEL NOTE] Note on the beveling process:

	Each of the bevel_vertex, bevel_edge, and extrude_face functions do not represent
	a full bevel/extrude operation. Instead, they should update the _connectivity_ of
	the mesh, _not_ the positions of newly created vertices. In fact, you should set
	the positions of new vertices to be exactly the same as wherever they "started from."

	When you click on a mesh element while in bevel mode, one of those three functions
	is called. But, because you may then adjust the distance/offset of the newly
	beveled face, we need another method of updating the positions of the new vertices.

	This is where bevel_positions and extrude_positions come in: these functions are
	called repeatedly as you move your mouse, the position of which determines the
	amount / shrink parameters. These functions are also passed an array of the original
	vertex positions, stored just after the bevel/extrude call, in order starting at
	face->halfedge->vertex, and the original element normal, computed just *before* the
	bevel/extrude call.

	Finally, note that the amount, extrude, and/or shrink parameters are not relative
	values -- you should compute a particular new position from them, not a delta to
	apply.
*/

/*
 * bevel_vertex: creates a face in place of a vertex
 *  v: the vertex to bevel
 *
 * returns: reference to the new face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_vertex(VertexRef v) {
	//A2Lx5 (OPTIONAL): Bevel Vertex
	// Reminder: This function does not update the vertex positions.
	// Remember to also fill in bevel_vertex_helper (A2Lx5h)

	(void)v;
    return std::nullopt;
}

/*
 * bevel_edge: creates a face in place of an edge
 *  e: the edge to bevel
 *
 * returns: reference to the new face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_edge(EdgeRef e) {
	//A2Lx6 (OPTIONAL): Bevel Edge
	// Reminder: This function does not update the vertex positions.
	// remember to also fill in bevel_edge_helper (A2Lx6h)

	(void)e;
    return std::nullopt;
}

/*
 * extrude_face: creates a face inset into a face
 *  f: the face to inset
 *
 * returns: reference to the inner face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::extrude_face(FaceRef f) {
    std::vector<HalfedgeRef> h;
    std::vector<VertexRef> v;
    std::vector<HalfedgeRef> t_orig;
    std::vector<EdgeRef> e_orig;

    HalfedgeRef cur = f->halfedge;
    do {
        h.push_back(cur);
        v.push_back(cur->vertex);
        t_orig.push_back(cur->twin);
        e_orig.push_back(cur->edge);
        cur = cur->next;
    } while (cur != f->halfedge);

    int n = (int)h.size();

    std::vector<VertexRef> vm(n);     
    std::vector<EdgeRef> e_inner(n);  
    std::vector<EdgeRef> e_spoke(n);  
    std::vector<FaceRef> quad(n);     
   
    std::vector<HalfedgeRef> a(n), b(n), c(n), d(n);

    for (int i = 0; i < n; i++) {
        vm[i] = emplace_vertex();
        vm[i]->position = v[i]->position;
        interpolate_data({v[i]}, vm[i]);

        e_inner[i] = emplace_edge();
        e_spoke[i] = emplace_edge();
        quad[i] = emplace_face(false);

        a[i] = emplace_halfedge();
        b[i] = emplace_halfedge();
        c[i] = emplace_halfedge();
        d[i] = emplace_halfedge();
    }

    for (int i = 0; i < n; i++) {
        int ni = (i + 1) % n;
        int pi = (i - 1 + n) % n;

        h[i]->vertex = vm[i];
        h[i]->twin = c[i];
        h[i]->edge = e_inner[i];
        h[i]->face = f;          
        e_inner[i]->halfedge = h[i];
        vm[i]->halfedge = h[i];

        a[i]->vertex = v[i];
        a[i]->twin = t_orig[i];
        a[i]->edge = e_orig[i];
        a[i]->face = quad[i];
        a[i]->next = b[i];
        t_orig[i]->twin = a[i];
        e_orig[i]->halfedge = a[i];

        b[i]->vertex = v[ni];
        b[i]->twin = d[ni];
        b[i]->edge = e_spoke[ni];
        b[i]->face = quad[i];
        b[i]->next = c[i];

        c[i]->vertex = vm[ni];
        c[i]->twin = h[i];
        c[i]->edge = e_inner[i];
        c[i]->face = quad[i];
        c[i]->next = d[i];

        d[i]->vertex = vm[i];
        d[i]->twin = b[pi];
        d[i]->edge = e_spoke[i];
        d[i]->face = quad[i];
        d[i]->next = a[i];
        e_spoke[i]->halfedge = d[i];

        quad[i]->halfedge = a[i];

        if (v[i]->halfedge == h[i]) {
            v[i]->halfedge = a[i];
        }
    }

    f->halfedge = h[0];
    return f;
}


/*
 * flip_edge: rotate non-boundary edge ccw inside its containing faces
 *  e: edge to flip
 *
 * if e is a boundary edge, does nothing and returns std::nullopt
 * if flipping e would create an invalid mesh, does nothing and returns std::nullopt
 *
 * otherwise returns the edge, post-rotation
 *
 * does not create or destroy mesh elements.
 */
std::optional<Halfedge_Mesh::EdgeRef> Halfedge_Mesh::flip_edge(EdgeRef e) {
    HalfedgeRef h = e->halfedge;
    HalfedgeRef t = h->twin;

    if (h->face->boundary || t->face->boundary) {
        return std::nullopt;
    }

    FaceRef F_h = h->face;
    FaceRef F_t = t->face;
    if (F_h == F_t) {
        return std::nullopt; 
    }

    VertexRef A = h->vertex;
    VertexRef B = t->vertex;
    HalfedgeRef hn = h->next; 
    HalfedgeRef tn = t->next; 
    VertexRef C = hn->next->vertex;
    VertexRef D = tn->next->vertex;

    if (C == D) {
        return std::nullopt; 
    }

    HalfedgeRef prev_h = hn;
    while (prev_h->next != h) {
        prev_h = prev_h->next;
    }

    HalfedgeRef prev_t = tn;
    while (prev_t->next != t) {
        prev_t = prev_t->next;
    }

    h->vertex = D;
    h->next   = hn->next;  
    t->vertex = C;
    t->next   = tn->next;  
    hn->next = t;
    hn->face = F_t;

    tn->next = h;
    tn->face = F_h;

    prev_h->next = tn; 
    prev_t->next = hn; 
    if (A->halfedge == h) {
        A->halfedge = tn;
    }
    if (B->halfedge == t) {
        B->halfedge = hn;
    }
    F_h->halfedge = h;
    F_t->halfedge = t;

    return e;
}


/*
 * make_boundary: add non-boundary face to boundary
 *  face: the face to make part of the boundary
 *
 * if face ends up adjacent to other boundary faces, merge them into face
 *
 * if resulting mesh would be invalid, does nothing and returns std::nullopt
 * otherwise returns face
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::make_boundary(FaceRef face) {
	//A2Lx7: (OPTIONAL) make_boundary

	return std::nullopt; //TODO: actually write this code!
}

/*
 * dissolve_vertex: merge non-boundary faces adjacent to vertex, removing vertex
 *  v: vertex to merge around
 *
 * if merging would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the merged face
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::dissolve_vertex(VertexRef v) {
    std::vector<HalfedgeRef> out_he;
    {
        HalfedgeRef cur = v->halfedge;
        do { out_he.push_back(cur); cur = cur->twin->next; } while (cur != v->halfedge);
    }
    int n = (int)out_he.size();

    for (auto h : out_he) {
        if (h->face->boundary) {
            return std::nullopt;
        }
    }

    std::vector<HalfedgeRef> end_he(n);
    for (int i = 0; i < n; i++) {
        HalfedgeRef t_prev = out_he[(i - 1 + n) % n]->twin;
        HalfedgeRef cur = out_he[i];
        while (cur->next != t_prev) cur = cur->next;
        end_he[i] = cur;
    }

    {
        int total = 0;
        for (int i = 0; i < n; i++) {
            HalfedgeRef cur = out_he[i]->next;
            HalfedgeRef stop = out_he[(i - 1 + n) % n]->twin;
            while (cur != stop) { total++; cur = cur->next; }
        }
        if (total < 3) {
            return std::nullopt;
        }
    }

    FaceRef merged = out_he[0]->face;

    for (int i = 0; i < n; i++) {
        VertexRef u_i = out_he[i]->twin->vertex;
        if (u_i->halfedge == out_he[i]->twin) {
            u_i->halfedge = end_he[i];
        }
    }

    for (int i = 0; i < n; i++) {
        end_he[i]->next = out_he[(i - 1 + n) % n]->next;
    }

    {
        HalfedgeRef start = out_he[0]->next;
        HalfedgeRef cur = start;
        do { 
            cur->face = merged; cur = cur->next; 
        } while (cur != start);
        merged->halfedge = start;
    }

    for (int i = 1; i < n; i++) {
        erase_face(out_he[i]->face);
    }

    for (int i = 0; i < n; i++) {
        EdgeRef ei = out_he[i]->edge;
        erase_halfedge(out_he[i]->twin);
        erase_halfedge(out_he[i]);
        erase_edge(ei);
    }
    erase_vertex(v);

    return merged;
}

/*
 * dissolve_edge: merge the two faces on either side of an edge
 *  e: the edge to dissolve
 *
 * merging a boundary and non-boundary face produces a boundary face.
 *
 * if the result of the merge would be an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the merged face.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::dissolve_edge(EdgeRef e) {
	// A2Lx2 (OPTIONAL): dissolve_edge

	//Reminder: use interpolate_data() to merge corner_uv / corner_normal data
	
    return std::nullopt;
}

/* collapse_edge: collapse edge to a vertex at its middle
 *  e: the edge to collapse
 *
 * if collapsing the edge would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the newly collapsed vertex
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_edge(EdgeRef e) {
    HalfedgeRef h = e->halfedge;
    HalfedgeRef t = h->twin;
    VertexRef A = h->vertex;  
    VertexRef B = t->vertex;  
    FaceRef F_h = h->face;
    FaceRef F_t = t->face;

    std::unordered_set<uint32_t> A_neighbors;
    {
        HalfedgeRef cur = A->halfedge;
        do {
            A_neighbors.insert(cur->twin->vertex->id);
            cur = cur->twin->next;
        } while (cur != A->halfedge);
    }
    std::unordered_set<uint32_t> face_h_verts, face_t_verts;
    if (!F_h->boundary) {
        HalfedgeRef cur = h->next;
        while (cur != h) { face_h_verts.insert(cur->vertex->id); cur = cur->next; }
    }
    if (!F_t->boundary) {
        HalfedgeRef cur = t->next;
        while (cur != t) { face_t_verts.insert(cur->vertex->id); cur = cur->next; }
    }
    {
        HalfedgeRef cur = B->halfedge;
        do {
            uint32_t nb = cur->twin->vertex->id;
            if (nb != A->id && A_neighbors.count(nb)) {
                if (!face_h_verts.count(nb) && !face_t_verts.count(nb))
                    return std::nullopt;
            }
            cur = cur->twin->next;
        } while (cur != B->halfedge);
    }

    auto count_sides = [](HalfedgeRef s) {
        int n = 0; HalfedgeRef c = s;
        do { n++; c = c->next; } while (c != s);
        return n;
    };
    bool h_tri = !F_h->boundary && count_sides(h) == 3;
    bool t_tri = !F_t->boundary && count_sides(t) == 3;

    HalfedgeRef h_next = h->next;
    HalfedgeRef t_next = t->next;
    HalfedgeRef h_prev = h_next; while (h_prev->next != h) h_prev = h_prev->next;
    HalfedgeRef t_prev = t_next; while (t_prev->next != t) t_prev = t_prev->next;

    VertexRef vm = emplace_vertex();
    vm->position = (A->position + B->position) / 2.0f;
    interpolate_data({A, B}, vm);

    {
        HalfedgeRef cur = t->next;
        while (cur != h) { cur->vertex = vm; cur = cur->twin->next; }
    }
    {
        HalfedgeRef cur = h->next;
        while (cur != t) { cur->vertex = vm; cur = cur->twin->next; }
    }

    if (h_tri) {
        HalfedgeRef hna = h_next;
        HalfedgeRef hnb = h_next->next;
        HalfedgeRef ta  = hna->twin;  
        HalfedgeRef tb  = hnb->twin;  

        EdgeRef hna_edge = hna->edge;
        EdgeRef hnb_edge = hnb->edge;

        ta->twin = tb;
        tb->twin = ta;
        ta->edge = hnb_edge;
        hnb_edge->halfedge = ta;

        VertexRef C = ta->vertex;  
        if (C->halfedge == hnb) {
            C->halfedge = ta;
        }

        if (ta->face->halfedge == hna) {
            ta->face->halfedge = ta;
        }
        if (tb->face->halfedge == hnb) {
            tb->face->halfedge = tb;
        }

        vm->halfedge = tb;

        erase_halfedge(hna);
        erase_edge(hna_edge);
        erase_halfedge(hnb);
        erase_face(F_h);
    } 
    else {
        h_prev->next = h_next;
        if (F_h->halfedge == h) {
            F_h->halfedge = h_next;
        }
        vm->halfedge = h_next;
    }

    if (t_tri) {
        HalfedgeRef tna = t_next;
        HalfedgeRef tnb = t_next->next;
        HalfedgeRef tc  = tna->twin;  
        HalfedgeRef td  = tnb->twin;  

        EdgeRef tna_edge = tna->edge;
        EdgeRef tnb_edge = tnb->edge;

        tc->twin = td;
        td->twin = tc;
        tc->edge = tnb_edge;
        tnb_edge->halfedge = tc;

        VertexRef D = tc->vertex;  
        if (D->halfedge == tnb) {
            D->halfedge = tc;
        }

        if (tc->face->halfedge == tna) {
            tc->face->halfedge = tc;
        }
        if (td->face->halfedge == tnb) {
            td->face->halfedge = td;
        }

        if (vm->halfedge == tna || vm->halfedge == tnb)
            vm->halfedge = td;

        erase_halfedge(tna);
        erase_edge(tna_edge);
        erase_halfedge(tnb);
        erase_face(F_t);
    } 
    else {
        t_prev->next = t_next;
        if (F_t->halfedge == t) {
            F_t->halfedge = t_next;
        }
        if (vm->halfedge == h) {
            vm->halfedge = t_next;
        }
    }

    if (A->halfedge == h) {
        A->halfedge = vm->halfedge;
    }
    if (B->halfedge == t) {
        B->halfedge = vm->halfedge;
    }

    erase_halfedge(h);
    erase_halfedge(t);
    erase_edge(e);
    erase_vertex(A);
    erase_vertex(B);

    return vm;
}

/*
 * collapse_face: collapse a face to a single vertex at its center
 *  f: the face to collapse
 *
 * if collapsing the face would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the newly collapsed vertex
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_face(FaceRef f) {
	//A2Lx3 (OPTIONAL): Collapse Face

	//Reminder: use interpolate_data() to merge corner_uv / corner_normal data on halfedges
	// (also works for bone_weights data on vertices!)

    return std::nullopt;
}

/*
 * weld_edges: glue two boundary edges together to make one non-boundary edge
 *  e, e2: the edges to weld
 *
 * if welding the edges would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns e, updated to represent the newly-welded edge
 */
std::optional<Halfedge_Mesh::EdgeRef> Halfedge_Mesh::weld_edges(EdgeRef e, EdgeRef e2) {
	//A2Lx8: Weld Edges

	//Reminder: use interpolate_data() to merge bone_weights data on vertices!

    return std::nullopt;
}



/*
 * bevel_positions: compute new positions for the vertices of a beveled vertex/edge
 *  face: the face that was created by the bevel operation
 *  start_positions: the starting positions of the vertices
 *     start_positions[i] is the starting position of face->halfedge(->next)^i
 *  direction: direction to bevel in (unit vector)
 *  distance: how far to bevel
 *
 * push each vertex from its starting position along its outgoing edge until it has
 *  moved distance `distance` in direction `direction`. If it runs out of edge to
 *  move along, you may choose to extrapolate, clamp the distance, or do something
 *  else reasonable.
 *
 * only changes vertex positions (no connectivity changes!)
 *
 * This is called repeatedly as the user interacts, just after bevel_vertex or bevel_edge.
 * (So you can assume the local topology is set up however your bevel_* functions do it.)
 *
 * see also [BEVEL NOTE] above.
 */
void Halfedge_Mesh::bevel_positions(FaceRef face, std::vector<Vec3> const &start_positions, Vec3 direction, float distance) {
	//A2Lx5h / A2Lx6h (OPTIONAL): Bevel Positions Helper
	
	// The basic strategy here is to loop over the list of outgoing halfedges,
	// and use the preceding and next vertex position from the original mesh
	// (in the start_positions array) to compute an new vertex position.
	
}

/*
 * extrude_positions: compute new positions for the vertices of an extruded face
 *  face: the face that was created by the extrude operation
 *  move: how much to translate the face
 *  shrink: amount to linearly interpolate vertices in the face toward the face's centroid
 *    shrink of zero leaves the face where it is
 *    positive shrink makes the face smaller (at shrink of 1, face is a point)
 *    negative shrink makes the face larger
 *
 * only changes vertex positions (no connectivity changes!)
 *
 * This is called repeatedly as the user interacts, just after extrude_face.
 * (So you can assume the local topology is set up however your extrude_face function does it.)
 *
 * Using extrude face in the GUI will assume a shrink of 0 to only extrude the selected face
 * Using bevel face in the GUI will allow you to shrink and increase the size of the selected face
 * 
 * see also [BEVEL NOTE] above.
 */
void Halfedge_Mesh::extrude_positions(FaceRef face, Vec3 move, float shrink) {
    std::vector<HalfedgeRef> inner_he;
    std::vector<Vec3> outer_pos;

    HalfedgeRef h = face->halfedge;
    do {
        inner_he.push_back(h);
        Vec3 op = h->twin->next->twin->vertex->position;
        outer_pos.push_back(op);
        h = h->next;
    } while (h != face->halfedge);

    int n = (int)inner_he.size();

    Vec3 centroid(0.0f, 0.0f, 0.0f);
    for (auto& p : outer_pos) {
        centroid += p;
    }
    centroid /= (float)n;

    for (int i = 0; i < n; i++) {
        VertexRef vm = inner_he[i]->vertex;
        vm->position = centroid + (outer_pos[i] - centroid) * (1.0f - shrink) + move;
    }
}
