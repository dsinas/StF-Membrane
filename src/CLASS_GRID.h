/*******************************************************************************
* 
* CLASS_GRID.h: A class for building a 2D lattice for membrane composition
* and the one for membrane out-of-plane height.
*
* This program is made freely available with the understanding that every 
* copy of this file must include this header and that it comes without any 
* WITHOUT ANY WARRANTY.
*
*******************************************************************************/

#ifndef CLASS_GRID_H
#define CLASS_GRID_H

using namespace std;

class LatticeClass {

protected:
    int m_lx, m_ly;

    /** declaring the spin lattice **/
    vector<int> m_spin;

    /** declaring the height lattice **/
    vector<double> m_height;

    /** nearest neighbors index matrices **/
    vector<vector<int> > m_neighbors;
    vector<vector<int> > m_next1_neighbors;
    vector<vector<int> > m_next2_neighbors;
    vector<vector<int> > m_next3_neighbors;
    vector<vector<int> > m_next4_neighbors;

public:

    LatticeClass(const int Lx, const int Ly): m_lx(Lx), m_ly(Ly),
            m_spin(Lx*Ly), m_height(Lx*Ly)
    {
        construct_lattices(m_lx,m_ly);
    }

    ~LatticeClass() {}

    /**
    convert the cell index (i,j) to 'idx' so that columns are filled one by one
    **/
    int get_cell_index(int i,int j) const {
        while (i<0) i += m_lx;
        while (i>=m_lx) i -= m_lx;
        while (j<0) j += m_ly;
        while (j>=m_ly) j -= m_ly;
        return i + j * m_lx;
    }

    /** convert the cell index idx to (i,j) (row & column respectively) **/
    vector<int> get_ij_index(const int idx) {
        if (idx>=get_volume() || idx<0) {
            cout<<"error: in get_ij_index() 'CLASS_GRID.h'"<<endl<<
                "index value exceeds the limits"<<endl;
            exit(1);
        }

        vector<int> result;
        result.resize(2);
        result[0] = idx % m_lx;
        result[1] = (int) (idx/m_lx);
        return result;
    }

    /** constructing the lattices and neighbors **/
    void construct_lattices(const int Lx,const int Ly) {
        m_spin.resize(Lx*Ly);
        m_height.resize(Lx*Ly);
        vector<int> fourSpins(4);
        vector<int> eightSpins(8);

        m_neighbors.resize(Lx*Ly,fourSpins);

        m_next1_neighbors.resize(Lx*Ly,fourSpins);
        m_next2_neighbors.resize(Lx*Ly,fourSpins);
        m_next3_neighbors.resize(Lx*Ly,eightSpins);
        m_next4_neighbors.resize(Lx*Ly,fourSpins);

        for (int i=0;i<Lx;i++)
            for (int j=0;j<Ly;j++) {
                int idx = get_cell_index(i,j);
                m_neighbors[idx].resize(4);
                m_neighbors[idx][0] = get_cell_index(i-1,j);
                m_neighbors[idx][1] = get_cell_index(i,j-1);
                m_neighbors[idx][2] = get_cell_index(i,j+1);
                m_neighbors[idx][3] = get_cell_index(i+1,j);
            }

        for (int i=0;i<Lx;i++)
            for (int j=0;j<Ly;j++) {
                int idx = get_cell_index(i,j);
                m_next1_neighbors[idx].resize(4);
                m_next1_neighbors[idx][0] = get_cell_index(i-1,j-1);
                m_next1_neighbors[idx][1] = get_cell_index(i-1,j+1);
                m_next1_neighbors[idx][2] = get_cell_index(i+1,j-1);
                m_next1_neighbors[idx][3] = get_cell_index(i+1,j+1);
            }

        for (int i=0;i<Lx;i++)
            for (int j=0;j<Ly;j++) {
                int idx = get_cell_index(i,j);
                m_next2_neighbors[idx].resize(4);
                m_next2_neighbors[idx][0] = get_cell_index(i-2,j);
                m_next2_neighbors[idx][1] = get_cell_index(i,j-2);
                m_next2_neighbors[idx][2] = get_cell_index(i,j+2);
                m_next2_neighbors[idx][3] = get_cell_index(i+2,j);
            }

        for (int i=0;i<Lx;i++)
            for (int j=0;j<Ly;j++) {
                int idx = get_cell_index(i,j);
                m_next3_neighbors[idx].resize(8);
                m_next3_neighbors[idx][0] = get_cell_index(i-2,j-1);
                m_next3_neighbors[idx][1] = get_cell_index(i-2,j+1);
                m_next3_neighbors[idx][2] = get_cell_index(i-1,j-2);
                m_next3_neighbors[idx][3] = get_cell_index(i-1,j+2);
                m_next3_neighbors[idx][4] = get_cell_index(i+1,j-2);
                m_next3_neighbors[idx][5] = get_cell_index(i+1,j+2);
                m_next3_neighbors[idx][6] = get_cell_index(i+2,j-1);
                m_next3_neighbors[idx][7] = get_cell_index(i+2,j+1);
            }

        for (int i=0;i<Lx;i++)
            for (int j=0;j<Ly;j++) {
                int idx = get_cell_index(i,j);
                m_next4_neighbors[idx].resize(4);
                m_next4_neighbors[idx][0] = get_cell_index(i-2,j-2);
                m_next4_neighbors[idx][1] = get_cell_index(i-2,j+2);
                m_next4_neighbors[idx][2] = get_cell_index(i+2,j-2);
                m_next4_neighbors[idx][3] = get_cell_index(i+2,j+2);
            }
    }

    inline void set_cellSpin(const int idx, const int _spin) {
        m_spin[idx] = _spin;
    }

    inline void flip_spin(const int idx) {
        m_spin[idx] *= -1;
    }

    inline void swap_spins(const int i, const int j) {
        int temp_spin = m_spin[i];
        m_spin[i] = m_spin[j];
        m_spin[j] = temp_spin;
    }

    inline int get_cellSpin(const int idx) {
        return m_spin[idx];
    }

    inline int get_cellSpin(const int i, const int j) {
        return get_cellSpin(get_cell_index(i,j));
    }

    inline int get_neighborsIndex(const int cell_index, const int neighbor) {
        return m_neighbors[cell_index][neighbor];
    }

    /** ******************************************************************* **/

    inline int get_sum_neighbors_spin(const int idx) const {
        const int* ptr = &(m_neighbors[idx][0]);
        return m_spin[ptr[0]] + m_spin[ptr[1]] + m_spin[ptr[2]] + m_spin[ptr[3]];
    }

    inline int get_sum_next1_neighbors_spin(const int idx) const {
        const int* ptr = &(m_next1_neighbors[idx][0]);
        return m_spin[ptr[0]] + m_spin[ptr[1]] + m_spin[ptr[2]] + m_spin[ptr[3]];
    }

    /** ******************************************************************* **/

    inline double get_sum_neighbors_height(const int idx) const {
        const int* ptr = &(m_neighbors[idx][0]);
        return m_height[ptr[0]] + m_height[ptr[1]] + m_height[ptr[2]] + m_height[ptr[3]];
    }

    inline double get_sum_next1_neighbors_height(const int idx) const {
        const int* ptr = &(m_next1_neighbors[idx][0]);
        return m_height[ptr[0]] + m_height[ptr[1]] + m_height[ptr[2]] + m_height[ptr[3]];
    }

    inline double get_sum_next2_neighbors_height(const int idx) const {
        const int* ptr = &(m_next2_neighbors[idx][0]);
        return m_height[ptr[0]] + m_height[ptr[1]] + m_height[ptr[2]] + m_height[ptr[3]];
    }

    inline double get_sum_next3_neighbors_height(const int idx) const {
        const int* ptr = &(m_next3_neighbors[idx][0]);
        return m_height[ptr[0]] + m_height[ptr[1]] + m_height[ptr[2]] + m_height[ptr[3]] +
               m_height[ptr[4]] + m_height[ptr[5]] + m_height[ptr[6]] + m_height[ptr[7]];
    }

    inline double get_sum_next4_neighbors_height(const int idx) const {
        const int* ptr = &(m_next4_neighbors[idx][0]);
        return m_height[ptr[0]] + m_height[ptr[1]] + m_height[ptr[2]] + m_height[ptr[3]];
    }

    /** get squared gradient at point (i,j) **/
    inline double get_delH2(const int idx) const {
        const int* ptr = &(m_neighbors[idx][0]);
        return 0.25*((m_height[ptr[0]]-m_height[ptr[3]])*(m_height[ptr[0]]-m_height[ptr[3]]) +
                     (m_height[ptr[1]]-m_height[ptr[2]])*(m_height[ptr[1]]-m_height[ptr[2]]));
    }

    /** ******************************************************************* **/
    /** declaring functions related to height of membrane ***************** **/

    inline void set_cellHeight(const int i, const int j, const double _height) {
        m_height[get_cell_index(i,j)] = _height;
    }

    inline void set_cellHeight(const int idx, const double _height) {
        m_height[idx] = _height;
    }

    inline double get_cellHeight(const int idx) {
        return m_height[idx];
    }

    inline double get_cellHeight(const int i, const int j) {
        return get_cellHeight(get_cell_index(i,j));
    }

    /** ******************************************************************* **/
    /** set & get system sizes ******************************************** **/

    void set_latticeSize(const int Lx, const int Ly) {
        m_lx = Lx;
        m_ly = Ly;
    }

    inline int get_lx() const {
        return m_lx;
    }

    inline int get_ly() const {
        return m_ly;
    }

    inline int get_volume() const {
        return (int) m_spin.size();
    }

    inline int get_neighbor_size() const {
        return (int) m_neighbors.at(0).size();
    }

};

#endif
