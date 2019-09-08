#ifndef __PLOT3D_H__
#define __PLOT3D_H__

#include <vector>

class PLOT3D_BLK
{
private:
    size_t m_nI, m_nJ, m_nK;
	bool m_3d;
	int m_dim;
	double *m_coord;
	size_t m_nIJ, m_nIJK;

public:
	PLOT3D_BLK(size_t nI, size_t nJ):
		m_nI(nI),
		m_nJ(nJ),
		m_nK(1),
		m_3d(false),
		m_dim(2),
		m_coord(new double[2*nI*nJ])
	{
		if(nI==0 || nJ == 0)
			throw("Invalid size.");

		m_nIJK = m_nIJ = m_nI * m_nJ;
	}
	
	PLOT3D_BLK(size_t nI, size_t nJ, size_t nK):
		m_nI(nI),
		m_nJ(nJ),
		m_nK(nK),
		m_3d(true),
		m_dim(3),
		m_coord(new double[3*nI*nJ*nK])
	{
		if(nI < 2 || nJ < 2 || nK < 2)
			throw("Invalid size.");

		m_nIJ = m_nI * m_nJ;
		m_nIJK = m_nI * m_nJ * m_nK;
	}

	~PLOT3D_BLK()
	{
		delete [] m_coord;
	}
	
	bool is3D() const
	{
		return m_3d;
	}

	int dimension() const
	{
		return m_dim;
	}

	size_t node_num() const
	{
		return m_nIJK;
	}

	size_t cell_num() const
	{
		size_t ret = 0;

		if(is3D)
			ret = (m_nI-1)*(m_nJ-1)*(m_nK-1);
		else
			ret = (m_nI-1)*(m_nJ-1);

		return ret;
	}

	size_t face_num() const
	{
		return boundary_face_num() + internal_face_num();
	}

	size_t boundary_face_num() const
	{
		size_t ret = 0;

		// TODO

		return ret;
	}

	size_t internal_face_num() const
	{
		size_t ret = 0;

		// TODO

		return ret;
	}

	size_t IMAX() const
	{
		return m_nI;
	}

	size_t JMAX() const
	{
		return m_nJ;
	}
	
	size_t KMAX() const
	{
		return m_nK;
	}

	// 0-based
	double *at(size_t i, size_t j) 
	{
		return m_coord + STX(i, j);
	}

	double *at(size_t i, size_t j, size_t k) 
	{
		return m_coord + STX(i, j, k);
	}

	// 1-based
	double *operator()(size_t i, size_t j) 
	{
		return at(i-1, j-1);
	}

	double *operator()(size_t i, size_t j, size_t k) 
	{
		return at(i-1, j-1, k-1);
	}

private:
	size_t STX(size_t i, size_t j) const
	{
		return i + j * m_nI;
	}

	size_t STX(size_t i, size_t j, size_t k) const
	{
		return i + j * m_nI + k * m_nIJ;
	}
};


class PLOT3D
{
private:
    size_t m_nBLK;
    bool m_3d;
	int m_dim;
	std::vector<PLOT3D_BLK*> m_blk;

public:
    PLOT3D() {}

	~PLOT3D() {}

	bool is3D() const
	{
		return m_3d;
	}

	int dimension() const
	{
		return m_dim;
	}

	size_t blk_num() const
	{
		return m_nBLK;
	}
};

#endif
