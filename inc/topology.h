#ifndef __TOPOLOGY_H__
#define __TOPOLOGY_H__

#include <cstddef>
#include <vector>

class Face;
class Cell;

class Node
{
private:
	size_t m_idx;
	double m_location[3];
	std::vector<size_t> m_nodeAdjoinIdx;
	std::vector<Node *> m_nodeAdjoinPtr;
	std::vector<size_t> m_faceBelongIdx;
	std::vector<Face *> m_faceBelongPtr;
	std::vector<size_t> m_cellBelongIdx;
	std::vector<Cell *> m_cellBelongPtr;

public:
	Node() :
		m_idx(0),
		m_location{ 0.0, 0.0, 0.0 }
	{
	}

	~Node() = default;

	double x() const
	{
		return m_location[0];
	}

	double &x()
	{
		return m_location[0];
	}

	double y() const
	{
		return m_location[1];
	}

	double &y()
	{
		return m_location[1];
	}

	double z() const
	{
		return m_location[2];
	}

	double &z()
	{
		return m_location[2];
	}
};

class Face 
{	
private:
	size_t m_idx;
	size_t m_cl, m_cr;
	std::vector<size_t> m_nodeWithinIdx;
	std::vector<Node *> m_nodeWithinPtr;
	double m_area;
	double m_center[3];
	double m_normalUnitVecLR[3], m_normalUnitVecRL[3];

public:
	Face()
	{
	}

	~Face() = default;

};

class Cell
{
private:
	size_t m_idx;
	double m_center[3];
	double m_volum;
	std::vector<size_t> m_nodeWithinIdx;
	std::vector<Node *> m_nodeWithinPtr;
	std::vector<size_t> m_faceWithinIdx;
	std::vector<Face *> m_faceWithinPtr;
	std::vector<size_t> m_adjCelIdx;
	std::vector<Cell *> m_adjCellnPtr;

public:
	Cell()
	{
	}

	~Cell() = default;
};

#endif
