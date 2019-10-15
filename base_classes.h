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
	bool GenMe();
	void MakeBasis();

	double bond_length, angle, dihedral;
	bool genme, rotatory; //genme is universally true except for the first methyl group
	int rotatory_index;
        int rotations;
	Eigen::Vector3d plane, to_define_plane;
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

class SM : public sp3
{
public:
	SM(Atom * parent, double dihedral = 0);
};

class C : public sp2
{
public:
	C(Atom * parent, double dihedral = 0);
};
class CR : public sp2
{
public:
	CR(Atom * parent, double dihedral = 0);
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
	Residue(){side_rotatory_bonds_2 = 0;};
	virtual Residue * Spawn(){return this;};
	virtual Atom* GenAll(Atom * junction){return new Atom;};
	virtual void NameMe(){};
	void Generate(int rotations);
	
	Atom * junction;
	std::vector<Atom> atoms;
	std::vector<Atom*> main_rotatory_atoms, side_rotatory_atoms_1, side_rotatory_atoms_2;
	int main_rotatory_bonds, side_rotatory_bonds_1, side_rotatory_bonds_2;
	std::string resname;
	std::string res1;
	int resid;
};
