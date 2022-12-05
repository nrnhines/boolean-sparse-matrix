NEURON { SUFFIX nothing }

VERBATIM

#include <vector>
#include <set>
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


using rank_t = std::size_t;  // can save space with std::uint32_t

// BSM is BooleanSparseMatrix
class BSM {
public:
    BSM(rank_t nrow, rank_t ncol);
    ~BSM();
    void setelm(rank_t irow, rank_t jcol, bool value);
    bool getelm(rank_t irow, rank_t jcol);
    rank_t nelm(); // return number of true elements
    rank_t nrow() { return nrow_; }
    rank_t ncol() { return ncol_; }
    void pr(); // print full matrix (should be small) for testing
    BSM* mul(BSM* that); // return this * that
    void add(BSM* that); // in place this + that
    Object* obj_;

private:
    BSM* format(bool); // return copy of matrix in row (true) or col (false) format
    std::vector<std::set<rank_t> >m_;
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

static Object** bsm_mul(void* v) {
    BSM* m1 = (BSM*) v;
    BSM* m2 = bsm_arg(1);
    BSM* m3 = m1->mul(m2);
    return bsm_temp_objvar(m3);
}

static Object** bsm_add(void* v) {
    BSM* m1 = (BSM*) v;
    BSM* m2 = bsm_arg(1);
    m1->add(m2);
    return bsm_temp_objvar(m1);
}

// all methods that should return double
static Member_func BSM_member_func[] {
    {"getelm", bsm_getelm},
    {"nelm", bsm_nelm},
    {"nrow", bsm_nrow},
    {"ncol", bsm_ncol},
    {0, 0}
};

// all methods that should return Object**
static Member_ret_obj_func BSM_obj_func[] {
    {"setelm", bsm_setelm},
    {"pr", bsm_pr},
    {"mul", bsm_mul},
    {"add", bsm_add},
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

static bool intersects(std::set<rank_t>& a, std::set<rank_t>& b) {
    // this relys on ordering
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

    return false;
}

BSM* BSM::mul(BSM* that) { // m_ik = this_ij * that_jk
    assert(ncol() == that->nrow())
    assert(nrow())
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

void BSM::add(BSM* that) { // inplace this + that
}

ENDVERBATIM
