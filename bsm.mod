: BooleanSparseMatrix(nrow, ncol)
: becomes available in Hoc (and Python) when
: register_BooleanSparseMatrix() is called.

NEURON { SUFFIX nothing }

VERBATIM

#include <vector>
#include <unordered_set>
// Perhaps the provision below for implementation with  ordered #include <set>
// can be removed since a timing test illustrated that
/**
unordered int32_t
BooleanSparseMatrix[0](10000, 10000) 99874 true elements 
BooleanSparseMatrix[1](10000, 10000) 100318 true elements
BooleanSparseMatrix[2](10000, 10000) 996647 true elements
matmul time 52.8998

ordered size_t
matmul time 54.0312

unordered size_t
matmul time 52.3568
**/
// It was originally thought that matmul would benefit from ordered sets.
// See bool intersected(...) below.

#if 0 // until classreg.h is installed
#include "classreg.h"
#else
extern void class2oc(const char*,
                     void* (*cons)(Object*),
                     void (*destruct)(void*),
                     Member_func*,
                     int (*checkpoint)(void**),
                     Member_ret_obj_func*,
                     Member_ret_str_func*);
#endif


using rank_t = std::size_t;  // if 4B enough, can save space with std::uint32_t
static bool unordered{true}; // but could easily change to ordered set
using rowset = std::unordered_set<rank_t>;

// BSM is BooleanSparseMatrix
class BSM {
public:
    BSM(rank_t nrow, rank_t ncol);
    ~BSM();
    void setelm(rank_t irow, rank_t jcol, bool value);
    bool getelm(rank_t irow, rank_t jcol);
    size_t nelm(); // return number of true elements
    rank_t nrow() { return nrow_; }
    rank_t ncol() { return ncol_; }
    void setrow(rank_t irow, IvocVect* indices);
    IvocVect* getrow(rank_t irow); // return HOC Vector of indices for that row.
    void pr(); // print full matrix (should be small) for testing
    BSM* clone();
    bool equal(BSM* that); // true if every element the same
    BSM* matmul(BSM* that); // return this * that
    void elmmul(BSM* that); // in place element wise this_ij * that_ij
    void add(BSM* that); // in place element wise this_ij + that_ij
    Object* obj_;

private:
    BSM* format(bool); // return copy of matrix in row (true) or col (false) format
    std::vector<rowset> m_;
    bool row_format_; // if true the sets are column indices.
    rank_t nrow_, ncol_;
};

static void* BSM_constructor(Object* obj) {
    rank_t nrow = (rank_t) chkarg(1, 0, 1e12);
    rank_t ncol = (rank_t) chkarg(2, 0, 1e12);
    BSM* m = new BSM(nrow, ncol);
    m->obj_ = obj;
    return (void*)m;
}

static void BSM_destructor(void* v) {
    delete (BSM*)v;
}

static Symbol* bsm_sym;

static BSM* bsm_arg(int i) {
    Object* ob = *hoc_objgetarg(i);
    if (!ob || ob->ctemplate != bsm_sym->u.ctemplate) {
        check_obj_type(ob, "BooleanSparseMatrix");
    }
    return (BSM*) (ob->u.this_pointer);
}

Object** bsm_temp_objvar(BSM* m) {
    Object** po;
    if (m->obj_) {
        po = hoc_temp_objptr(m->obj_);
    } else {
        po = hoc_temp_objvar(bsm_sym, (void*) m);
        m->obj_ = *po;   
    }
    return po;
}

static double bsm_getelm(void* v) {
    BSM* m = (BSM*) v;
    rank_t irow = (rank_t) chkarg(1, 0, m->nrow() - 1);
    rank_t jcol = (rank_t) chkarg(2, 0, m->ncol() - 1);
    return (double) m->getelm(irow, jcol);
}

static double bsm_nelm(void* v) {
    BSM* m = (BSM*) v;
    return (double) m->nelm();
}

static double bsm_nrow(void* v) {
    BSM* m = (BSM*) v;
    return (double) m->nrow();
}

static double bsm_ncol(void* v) {
    BSM* m = (BSM*) v;
    return (double) m->ncol();
}

static double bsm_equal(void* v) {
    BSM* m = (BSM*) v;
    BSM* m2 = bsm_arg(1);
    return (double) m->equal(m2);
}

static Object** bsm_setelm(void* v) {
    BSM* m = (BSM*) v;
    rank_t irow = (rank_t) chkarg(1, 0, m->nrow() - 1);
    rank_t jcol = (rank_t) chkarg(2, 0, m->ncol() - 1);
    bool val = (bool) chkarg(3, 0, 1);    
    m->setelm(irow, jcol, val);
    return bsm_temp_objvar(m);
}

static Object** bsm_pr(void* v) {
    BSM* m = (BSM*) v;
    m->pr();
    return bsm_temp_objvar(m);
}

static Object** bsm_clone(void* v) {
    BSM* m1 = (BSM*) v;
    return bsm_temp_objvar(m1->clone());
}

static Object** bsm_matmul(void* v) {
    BSM* m1 = (BSM*) v;
    BSM* m2 = bsm_arg(1);
    BSM* m3 = m1->matmul(m2);
    return bsm_temp_objvar(m3);
}

static Object** bsm_elmmul(void* v) {
    BSM* m1 = (BSM*) v;
    BSM* m2 = bsm_arg(1);
    m1->elmmul(m2);
    return bsm_temp_objvar(m1);
}

static Object** bsm_add(void* v) {
    BSM* m1 = (BSM*) v;
    BSM* m2 = bsm_arg(1);
    m1->add(m2);
    return bsm_temp_objvar(m1);
}

static Object** bsm_getrow(void* v) {
    BSM* m = (BSM*)v;
    rank_t irow = (rank_t) chkarg(1, 0, m->nrow() - 1);
    IvocVect* vrow = m->getrow(irow);
    return vector_temp_objvar(vrow);
}

static Object** bsm_setrow(void* v) {
    BSM* m = (BSM*)v;
    rank_t irow = (rank_t) chkarg(1, 0, m->nrow() - 1);
    IvocVect* vrow = vector_arg(2);
    m->setrow(irow, vrow);
    return bsm_temp_objvar(m);
}

// all methods that should return double
static Member_func BSM_member_func[] {
    {"getelm", bsm_getelm},
    {"nelm", bsm_nelm},
    {"nrow", bsm_nrow},
    {"ncol", bsm_ncol},
    {"equal", bsm_equal},
    {0, 0}
};

// all methods that should return Object**
static Member_ret_obj_func BSM_obj_func[] {
    {"setelm", bsm_setelm},
    {"pr", bsm_pr},
    {"clone", bsm_clone},
    {"matmul", bsm_matmul},
    {"elmmul", bsm_elmmul},
    {"add", bsm_add},
    {"getrow", bsm_getrow},
    {"setrow", bsm_setrow},
    {0, 0}
};

static void BSM_reg() {
    // once only
    static bool registered{false};
    if (registered) {
        return;
    }
    class2oc("BooleanSparseMatrix", BSM_constructor, BSM_destructor,
        BSM_member_func, nullptr, BSM_obj_func, nullptr);
    bsm_sym = hoc_lookup("BooleanSparseMatrix");
}
ENDVERBATIM

PROCEDURE register_BooleanSparseMatrix() {
VERBATIM
    BSM_reg();
ENDVERBATIM
}

VERBATIM

// The rest is c++ implementaton of BSM

BSM::BSM(rank_t nrow, rank_t ncol) {
    nrow_ = nrow;
    ncol_ = ncol;
    row_format_ = true;
    m_.resize(nrow_);
    obj_ = nullptr;
}

BSM::~BSM() {
}

void BSM::setelm(rank_t irow, rank_t jcol, bool value) {
    if (value) {
        m_[irow].insert(jcol);
    } else {
        m_[irow].erase(jcol);
    }
}

bool BSM::getelm(rank_t irow, rank_t jcol) {
    if (m_[irow].count(jcol)) {
        return true;
    }
    return false;
}

std::size_t BSM::nelm() {
    std::size_t n{0};
    for(auto& s: m_) {
        n += s.size();
    }
    return n;
}

bool BSM::equal(BSM* m2) {
    if (nrow() != m2->nrow() || ncol() != m2->ncol()) {
        return false;
    }
    for (rank_t i = 0; i < nrow(); ++i) {
        if (m_[i] != m2->m_[i]) {
            return false;
        }
    }
    return true;
}

void BSM::pr() {
    Printf("%s(%zd, %zd) %zd true elements\n", hoc_object_name(obj_), (size_t)nrow_, (size_t)ncol_, nelm());
    auto nr = (nrow_ < 50) ? nrow_ : 50;
    auto nc = (ncol_ < 50) ? ncol_ : 50;
    for (rank_t i = 0; i < nr; ++i) {
        for (rank_t j = 0; j < nc; ++j) {
            Printf("%s", getelm(i, j) ? "X" : ".");
        }
        Printf("\n");
    }
}

static bool intersects(rowset& a, rowset& b) {
    if (unordered) {
        for(auto ia: a) { // may want to interate over set with fewer elements
            if (b.count(ia)) {
                return true;
            }
        }
    } else {
        // this relies on ordering
        auto ia = a.begin();
        auto ib = b.begin();

        while(ia != a.end() && ib != b.end()) {
            if (*ia == *ib) {
                return true;
            }
            if (*ia < *ib) {
                ++ia;
            } else {
                ++ib;
            }
        }
    }

    return false;
}

BSM* BSM::matmul(BSM* that) { // m_ik = this_ij * that_jk
    assert(ncol() == that->nrow())
    BSM* m2 = that->format(false); // column format
    BSM* m = new BSM(nrow_, that->ncol());
    for (rank_t i = 0; i < nrow(); ++i) {
        auto& this_row_set = m_[i];
        for (rank_t k = 0; k < that->ncol(); ++k) {
            auto& that_col_set = m2->m_[k];
            m->setelm(i, k, intersects(this_row_set, that_col_set));
        }
    }
    delete m2;
    return m;
}

// return copy of matrix in row (true) or col (false) format
BSM* BSM::format(bool format) {
    assert(row_format_);
    BSM* m = new BSM(nrow(), ncol());
    m->row_format_ = false;
    m->m_.resize(ncol());
    for (rank_t irow = 0; irow < nrow(); ++irow) {
        auto& s = m_[irow];
        for (auto& jcol: s) {
            m->m_[jcol].insert(irow);
        }
    }
    return m;
}

BSM* BSM::clone() { // deep copy
    BSM* m = new BSM(nrow(), ncol());
    m->m_ = m_;
    return m;
}

void BSM::elmmul(BSM* that) { // inplace element wise this_ij * that_ij
    assert(nrow() == that->nrow());
    assert(ncol() == that->ncol());
    for (rank_t irow = 0; irow < nrow(); ++irow) {
        auto& row = m_[irow];
        auto row1 = row; // copy because will change row (really needed?)
        auto& row2 = that->m_[irow];
        for (auto j: row1) {
            if (row2.count(j) == 0) {
                row.erase(j);
            }
        }
    }
}

void BSM::add(BSM* that) { // inplace this + that
    assert(nrow() == that->nrow());
    assert(ncol() == that->ncol());
    for (rank_t irow = 0; irow < nrow(); ++irow) {
        auto& row = m_[irow];
        auto& row2 = that->m_[irow];
        for (auto j: row2) {
            row.insert(j);
        }
    }

}

IvocVect* BSM::getrow(rank_t irow) {
    assert(irow < nrow());
    auto& row = m_[irow];
    IvocVect* vrow = vector_new1(row.size());
    double* px;
    vector_instance_px(vrow, &px);
    size_t iv{0};
    for (auto j: row) {
        px[iv++] = double(j);
    }
    return vrow;
}

void BSM::setrow(rank_t irow, IvocVect* vrow) {
    assert(irow < nrow());
    double* px;
    size_t nv = vector_instance_px(vrow, &px);
    auto& row = m_[irow];
    row.clear();
    for (size_t iv = 0; iv < nv; ++iv) {
        rank_t j = (rank_t)px[iv];
        assert(j < ncol());
        row.insert(j);
    }
}

ENDVERBATIM
