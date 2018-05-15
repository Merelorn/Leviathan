#pragma once
#include <eigen/Eigen/Geometry>
#include <vector>
#include <fstream>

class Atom
{
public:
	Atom();
	Atom(Eigen::Vector3d coords);
	Atom(Atom * parent, double angle, double dihedral = 0);
	Atom(Atom * parent, double dihedral);

	virtual void operator= (const Atom & rhs);
	void PrintMe(std::ofstream & myofile, int atomid, std::string resname, int resid);
	bool GenMe(int rotations);
	void MakeBasis(Eigen::Vector3d plane, Eigen::Vector3d parent_in_plane, Eigen::Vector3d & b1, Eigen::Vector3d & b2, double rotation, int rotations);

	double bond_length, angle, dihedral;
	bool genme, rotatory; //genme is universally true except for the first methyl group
	int rotatory_index;
	std::string type, pdb_type;
	Atom * parent, * plane_def;
	Eigen::Vector3d coordinates, b1, b2;
	std::vector<Atom*> children;
};

/// hybridizations

class sp3 : public Atom
{
public:
	sp3(Atom * parent, double dihedral = 0);
};

class sp2 : public Atom
{
public:
	sp2(Atom * parent, double dihedral = 0);
};

//// ATOM TYPES /////
class CA : public sp3
{
public:
	CA(Atom * parent, double dihedral = 0);
};

class H : public sp3
{
public:
	H(Atom * parent, double dihedral = 0);
};

class NZ : public sp3
{
public:
	NZ(Atom * parent, double dihedral = 0);
};

class SG : public sp3
{
public:
	SG(Atom * parent, double dihedral = 0);
};


class C : public sp2
{
public:
	C(Atom * parent, double dihedral = 0);
};
class O : public sp2
{
public:
	O(Atom * parent, double dihedral = 0);
};
class N : public sp2
{
public:
	N(Atom * parent, double dihedral = 0);
};
///////////////////


class Residue
{
public:
	Residue(){};
	virtual Residue * Spawn(){return this;};
	virtual Atom* GenAll(Atom * junction){return new Atom;};
	virtual void NameMe(){};
	void Generate(int rotations);
	
	Atom * junction;
	std::vector<Atom> atoms;
	std::vector<Atom*> rotatory_atoms;
	int rotatory_bonds;
	std::string resname;
	int resid;
};
