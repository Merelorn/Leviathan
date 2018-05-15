#include "base_classes.h"
#include <iostream>
#include <fstream>
#include <iomanip>

Atom::Atom(){parent = NULL; plane_def = NULL; rotatory = false;rotatory_index = 0; genme = true; genme = true; };
Atom::Atom(Eigen::Vector3d coords) : coordinates(coords){plane_def = NULL; parent = NULL;rotatory = false;rotatory_index = 0; genme = true; };
Atom::Atom(Atom * parent, double angle, double dihedral) : parent(parent), angle(angle), dihedral(dihedral){ plane_def = NULL; rotatory = false; rotatory_index = 0; genme = true; };
Atom::Atom(Atom * parent, double dihedral) : parent(parent), dihedral(dihedral){ plane_def = NULL; rotatory = false; rotatory_index = 0; genme = true; };
void Atom::operator=(const Atom & rhs)
{
	bond_length = rhs.bond_length;
	angle = rhs.angle;
	dihedral = rhs.dihedral;
	rotatory = rhs.rotatory;
	rotatory_index = rhs.rotatory_index;
	type = rhs.type;
	genme = rhs.genme;
	parent = rhs.parent;
	plane_def = rhs.plane_def;
	coordinates = rhs.coordinates;
	b1 = rhs.b1;
	b2 = rhs.b2;
	children = rhs.children;

	if ( parent != NULL )
	{
		parent->children.push_back(this);
	}
}


void Atom::PrintMe(std::ofstream & myofile, int atomid, std::string resname, int resid)
{
    myofile << "ATOM";
	myofile.width(7); myofile << atomid;
	myofile << "  ";
	myofile.width(3); myofile << pdb_type;
	myofile << " ";
	myofile.width(3); myofile << resname;
	myofile.width(6); myofile << resid;
	myofile << "    ";
	myofile << std::setprecision(3) << std::fixed;
	myofile.width(8); myofile << coordinates.x();
	myofile.width(8); myofile << coordinates.y();
	myofile.width(8); myofile << coordinates.z();
	myofile << "  1.00  0.00";
	myofile << "           ";
	myofile.width(1); myofile << type;
	myofile << std::endl;
}

bool Atom::GenMe(int rotations)
{
//	if ( b_angle && parent->parent == NULL ){return false;}
//	if ( b_dihedral && parent->parent->parent == NULL ){return false;}	// does not work?!!?!
//	if ( b_dihedral && ! b_angle ){return false;}

	double converted_dihedral = dihedral * 2 * M_PI / 360;

	double ccc = cos(parent->angle);

	coordinates = (parent->parent->coordinates - parent->coordinates);
	coordinates.normalize();
	coordinates *= ccc;
	coordinates += sqrt(1 - pow(ccc,2)) * (parent->b1 * cos(converted_dihedral) + parent->b2 * sin(converted_dihedral));
	coordinates *= bond_length;
	coordinates += parent->coordinates;

	Eigen::Vector3d plane, to_define_plane;
	plane = parent->coordinates - coordinates;
	to_define_plane = plane;

//	if ( parent->parent != NULL ){
		to_define_plane = parent->parent->coordinates - parent->coordinates;
/*		if ( parent->plane_def != NULL ){
			to_define_plane = parent->plane_def->coordinates - parent->coordinates;
		}
/*
	} else
*/
	MakeBasis(plane,to_define_plane,b1,b2,rotatory_index,rotations);

	return true;
}

void Atom::MakeBasis(Eigen::Vector3d plane, Eigen::Vector3d to_define_plane, Eigen::Vector3d & b1, Eigen::Vector3d & b2,double rotation, int rotations){

	Eigen::Vector3d bb1,bb2;

	bb2 = plane.cross(to_define_plane);
	if ( bb2.isZero() ){ bb2 = plane.cross(Eigen::Vector3d(1,0,0)); }
	if ( bb2.isZero() ){ bb2 = plane.cross(Eigen::Vector3d(0,1,0)); std::cout << "bump!" << std::endl; }		// !!!! not smooth !!!!
	// smoothness is not possible in odd-paired dimensions
	bb1 = plane.cross(bb2);
	
	bb1.normalize();
	bb2.normalize();

	rotation = 2 * M_PI * rotation / rotations;

	b1 = cos(rotation) * bb1 + sin(rotation) * bb2;
	b2 = -sin(rotation) * bb1 + cos(rotation) * bb2;

}


sp3::sp3(Atom * parent, double dihedral) : Atom(parent, (109.5 * M_PI / 180), dihedral){};

sp2::sp2(Atom * parent, double dihedral) : Atom(parent, (120 * M_PI / 180), dihedral){};


//// ATOM TYPES /////
CA::CA(Atom * parent, double dihedral) : sp3(parent, dihedral){bond_length = 1.5; type="C"; };
H::H(Atom * parent, double dihedral) : sp3(parent, dihedral){bond_length = 1.0; type="H"; };
C::C(Atom * parent, double dihedral) : sp2(parent, dihedral){bond_length = 1.5; type="C"; };
O::O(Atom * parent, double dihedral) : sp2(parent, dihedral){bond_length = 1.4; type="O"; };
N::N(Atom * parent, double dihedral) : sp2(parent, dihedral){bond_length = 1.4; type="N"; };
NZ::NZ(Atom * parent, double dihedral) : sp3(parent, dihedral){bond_length = 1.5; type="N"; };
SG::SG(Atom * parent, double dihedral) : sp3(parent, dihedral){bond_length = 1.8; type="S"; };
