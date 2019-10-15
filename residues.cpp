
#include "factory.h"

class Ccap0 : public Residue
{
public:

	Ccap0(){std::vector<Atom> dummy(6,Atom()); atoms = dummy;side_rotatory_bonds = 0; main_rotatory_bonds = 0;resname = "CAP";res1 = "";}

	Atom* GenAll(Atom * junction_atom)
	{
		bool test = true;
		Eigen::Vector3d c(0.000,0.000,0.000), h1(0.000,0.000,1.000);
		Atom cc(c);
		cc.angle = M_PI * 109.5 / 180;
		cc.type = "C";
		cc.rotations = 3;

		Atom hh(h1);
		hh.type = "H";
		hh.rotations = 3;

		atoms[0] = hh;  //atoms[0]
		atoms[0].plane = hh.coordinates - cc.coordinates;
		atoms[0].to_define_plane = hh.coordinates - cc.coordinates;
		atoms[0].MakeBasis();
		atoms[0].children.push_back(&atoms[1]);
		atoms[0].genme = false;

		cc.parent = &atoms[0];
		atoms[1] = cc;
		atoms[1].plane = hh.coordinates - cc.coordinates;
		atoms[1].to_define_plane = hh.coordinates - cc.coordinates;
		atoms[1].MakeBasis();
		atoms[1].genme = false;

		atoms[2] = H(&atoms[1],60.d);
		atoms[3] = H(&atoms[1],180.d);
		atoms[4] = C(&atoms[1],300.d); atoms[4].plane_def = &atoms[5];
		atoms[5] = O(&atoms[4],0.d);

		junction = &(*----atoms.end()); return junction; //penultimate member -> carbonyl C
	}
	void NameMe(){
		atoms[0].pdb_type = "1HA";
		atoms[1].pdb_type = "CA";
		atoms[2].pdb_type = "2HA";
		atoms[3].pdb_type = "3HA";
		atoms[4].pdb_type = "C";
		atoms[5].pdb_type = "O";
	}

};
static factory::Registrar<Ccap0> Ccap0_dummy("Ccap");	

class A2 : public Residue
{
public:

	A2(){std::vector<Atom> dummy(10,Atom()); atoms = dummy; side_rotatory_bonds = 0; main_rotatory_bonds = 2;resname = "ALA";res1 = "A";}
	Residue* Spawn(){return new A2;}

	Atom* GenAll(Atom * junction_atom)
	{
		atoms[0] = N(junction_atom,180.d);	// peptide nitrogen has dihedral 180, peptide oxygen has dihedral 0 (arbitrary choice)
		atoms[1] = H(&atoms[0],180.d);	// trans to oxygen
		atoms[2] = CA(&atoms[0],0.d);	main_rotatory_atoms.push_back(&atoms[2]); // trans
		atoms[3] = H(&atoms[2],180.d);

		atoms[4] = CA(&atoms[2],300.d);	//CB
		atoms[5] = H(&atoms[4],0.d);
		atoms[6] = H(&atoms[4],120.d);
		atoms[7] = H(&atoms[4],240.d);

		atoms[8] = C(&atoms[2],60.d); atoms[8].plane_def = &atoms[9]; main_rotatory_atoms.push_back(&atoms[8]); 
		atoms[9] = O(&atoms[8],0.d);

		junction = &(*----atoms.end()); return junction; //penultimate member -> carbonyl C
	}

	void NameMe(){
		atoms[0].pdb_type = "N";
		atoms[1].pdb_type = "H";
		atoms[2].pdb_type = "CA";
		atoms[3].pdb_type = "HA";
		atoms[4].pdb_type = "CB";
		atoms[5].pdb_type = "1HB";
		atoms[6].pdb_type = "2HB";
		atoms[7].pdb_type = "3HB";
		atoms[8].pdb_type = "C";
		atoms[9].pdb_type = "O";
	}

};
static factory::Registrar<A2> A2_dummy("A");

class C3 : public Residue
{
public:

	C3(){std::vector<Atom> dummy(11,Atom()); atoms = dummy;side_rotatory_bonds = 1;main_rotatory_bonds = 2;resname = "CYS";res1 = "C";}
	Residue* Spawn(){return new C3;}

		Atom* GenAll(Atom * junction_atom)
	{
		atoms[0] = N(junction_atom,180.d);	// peptide nitrogen has dihedral 180, peptide oxygen has dihedral 0 (arbitrary choice)
		atoms[1] = H(&atoms[0],180.d);	// trans to oxygen
		atoms[2] = CA(&atoms[0],0.d);	main_rotatory_atoms.push_back(&atoms[2]); // trans
		atoms[3] = H(&atoms[2],180.d);

		atoms[4] = CA(&atoms[2],300.d);	side_rotatory_atoms.push_back(&atoms[4]); //CB
		atoms[5] = H(&atoms[4],0.d);
		atoms[6] = H(&atoms[4],120.d);
		
		atoms[7] = SG(&atoms[4],240.d);	//SG
		atoms[8] = H(&atoms[7],0.d);

		atoms[9] = C(&atoms[2],60.d);	atoms[9].plane_def = &atoms[10]; main_rotatory_atoms.push_back(&atoms[9]); 
		atoms[10] = O(&atoms[9],0.d);

		junction = &(*----atoms.end()); return junction;
	}
	void NameMe(){
		atoms[0].pdb_type = "N";
		atoms[1].pdb_type = "H";
		atoms[2].pdb_type = "CA";
		atoms[3].pdb_type = "HA";
		atoms[4].pdb_type = "CB";
		atoms[5].pdb_type = "1HB";
		atoms[6].pdb_type = "2HB";
		atoms[7].pdb_type = "SG";
		atoms[8].pdb_type = "HG";
		atoms[9].pdb_type = "C";
		atoms[10].pdb_type = "O";
	}
};
static factory::Registrar<C3> C3_dummy("C");

class D4 : public Residue
{
public:

	D4(){std::vector<Atom> dummy(12,Atom()); atoms = dummy;main_rotatory_bonds = 2;side_rotatory_bonds=2;resname = "ASP";res1 = "D";}
	Residue* Spawn(){return new D4;}

	Atom* GenAll(Atom * junction_atom)
	{
		atoms[0] = N(junction_atom,180.d);	// peptide nitrogen has dihedral 180, peptide oxygen has dihedral 0 (arbitrary choice)
		atoms[1] = H(&atoms[0],180.d);	// trans to oxygen
		atoms[2] = CA(&atoms[0],0.d);	main_rotatory_atoms.push_back(&atoms[2]); // trans
		atoms[3] = H(&atoms[2],180.d);

		atoms[4] = CA(&atoms[2],300.d);	side_rotatory_atoms.push_back(&atoms[4]);	//CB
		atoms[5] = H(&atoms[4],0.d);
		atoms[6] = H(&atoms[4],120.d);

		atoms[7] = C(&atoms[4],240.d);	side_rotatory_atoms.push_back(&atoms[7]);	//CD
		atoms[8] = O(&atoms[7],0.d); atoms[7].plane_def = &atoms[8];
		atoms[9] = O(&atoms[7],180.d);

		atoms[10] = C(&atoms[2],60.d);	atoms[10].plane_def = &atoms[11];	main_rotatory_atoms.push_back(&atoms[10]);
		atoms[11] = O(&atoms[10],0.d);

		junction = &(*----atoms.end()); return junction;
	}
	void NameMe(){
		atoms[0].pdb_type = "N";
		atoms[1].pdb_type = "H";
		atoms[2].pdb_type = "CA";
		atoms[3].pdb_type = "HA";
		atoms[4].pdb_type = "CB";
		atoms[5].pdb_type = "1HB";
		atoms[6].pdb_type = "2HB";
		atoms[7].pdb_type = "CG";
		atoms[8].pdb_type = "OD1";
		atoms[9].pdb_type = "OD2";
		atoms[10].pdb_type = "C";
		atoms[11].pdb_type = "O";
	}
};
static factory::Registrar<D4> D4_dummy("D");

class E5 : public Residue
{
public:

	E5(){std::vector<Atom> dummy(15,Atom()); atoms = dummy;main_rotatory_bonds = 2;side_rotatory_bonds=3;resname = "GLU";res1 = "E";}
	Residue* Spawn(){return new E5;}

	Atom* GenAll(Atom * junction_atom)
	{
		atoms[0] = N(junction_atom,180.d);	// peptide nitrogen has dihedral 180, peptide oxygen has dihedral 0 (arbitrary choice)
		atoms[1] = H(&atoms[0],180.d);	// trans to oxygen
		atoms[2] = CA(&atoms[0],0.d);	main_rotatory_atoms.push_back(&atoms[2]); // trans
		atoms[3] = H(&atoms[2],180.d);

		atoms[4] = CA(&atoms[2],300.d);	side_rotatory_atoms.push_back(&atoms[4]);	//CB
		atoms[5] = H(&atoms[4],0.d);
		atoms[6] = H(&atoms[4],120.d);

		atoms[7] = CA(&atoms[4],240.d);	side_rotatory_atoms.push_back(&atoms[7]);	//CG
		atoms[8] = H(&atoms[7],0.d);
		atoms[9] = H(&atoms[7],120.d);

		atoms[10] = C(&atoms[7],240.d);	side_rotatory_atoms.push_back(&atoms[10]);	//CD
		atoms[11] = O(&atoms[10],0.d); atoms[10].plane_def = &atoms[11];
		atoms[12] = O(&atoms[10],180.d);

		atoms[13] = C(&atoms[2],60.d);	atoms[13].plane_def = &atoms[14];	main_rotatory_atoms.push_back(&atoms[13]);
		atoms[14] = O(&atoms[13],0.d);

		junction = &(*----atoms.end()); return junction;
	}
	void NameMe(){
		atoms[0].pdb_type = "N";
		atoms[1].pdb_type = "H";
		atoms[2].pdb_type = "CA";
		atoms[3].pdb_type = "HA";
		atoms[4].pdb_type = "CB";
		atoms[5].pdb_type = "1HB";
		atoms[6].pdb_type = "2HB";
		atoms[7].pdb_type = "CG";
		atoms[8].pdb_type = "1HG";
		atoms[9].pdb_type = "2HG";
		atoms[10].pdb_type = "CD";
		atoms[11].pdb_type = "OE1";
		atoms[12].pdb_type = "OE2";
		atoms[13].pdb_type = "C";
		atoms[14].pdb_type = "O";
	}
};
static factory::Registrar<E5> E5_dummy("E");

class F4 : public Residue
{
public:

	F4(){std::vector<Atom> dummy(20,Atom()); atoms = dummy;main_rotatory_bonds = 2;side_rotatory_bonds=2;resname = "PHE";res1 = "F";}
	Residue* Spawn(){return new F4;}

	Atom* GenAll(Atom * junction_atom)
	{
		atoms[0] = N(junction_atom,180.d);	// peptide nitrogen has dihedral 180, peptide oxygen has dihedral 0 (arbitrary choice)
		atoms[1] = H(&atoms[0],180.d);	// trans to oxygen
		atoms[2] = CA(&atoms[0],0.d);	main_rotatory_atoms.push_back(&atoms[2]); // trans
		atoms[3] = H(&atoms[2],180.d);

		atoms[4] = CA(&atoms[2],300.d);	side_rotatory_atoms.push_back(&atoms[4]);	//CB
		atoms[5] = H(&atoms[4],0.d);
		atoms[6] = H(&atoms[4],120.d);

		//aromatic_plane.coordinates = atoms[5]-atoms[4] + atoms[6]-atoms[4];

		atoms[7] = C(&atoms[4],240.d);	atoms[7].plane_def = atoms[7].parent;		side_rotatory_atoms.push_back(&atoms[7]);	//CG
		atoms[8] = C(&atoms[7],0.d);		atoms[8].plane_def = atoms[8].parent;		//CD
		atoms[9] = C(&atoms[7],180.d);	atoms[9].plane_def = atoms[9].parent;		//CD

		atoms[10] = C(&atoms[8],0.d);											//CE
		atoms[11] = C(&atoms[9],0.d);		atoms[11].plane_def = atoms[11].parent;	//CE

		atoms[12] = C(&atoms[11],180.d);										//CF

		atoms[12].children.push_back(&atoms[10]);
		atoms[10].children.push_back(&atoms[12]);

		atoms[13] = H(&atoms[8],180.d - atoms[10].dihedral);	//aromatic hydrogens
		atoms[14] = H(&atoms[9],180.d - atoms[11].dihedral);	//aromatic hydrogens
		atoms[15] = H(&atoms[10],0.d);	//aromatic hydrogens
		atoms[16] = H(&atoms[11],180.d - atoms[12].dihedral);	//aromatic hydrogens
		atoms[17] = H(&atoms[12],0.d);	//aromatic hydrogens

		
		
		atoms[18] = C(&atoms[2],60.d);	atoms[18].plane_def = &atoms[19];	main_rotatory_atoms.push_back(&atoms[18]);
		atoms[19] = O(&atoms[18],0.d);

		junction = &(*----atoms.end()); return junction;
	}
	void NameMe(){
		atoms[0].pdb_type = "N";
		atoms[1].pdb_type = "H";
		atoms[2].pdb_type = "CA";
		atoms[3].pdb_type = "HA";
		atoms[4].pdb_type = "CB";
		atoms[5].pdb_type = "1HB";
		atoms[6].pdb_type = "2HB";
		atoms[7].pdb_type = "CG";
		atoms[8].pdb_type = "CD1";
		atoms[9].pdb_type = "CD2";
		atoms[10].pdb_type = "CE1";
		atoms[11].pdb_type = "CE2";
		atoms[12].pdb_type = "CZ";
		atoms[13].pdb_type = "HD1";
		atoms[14].pdb_type = "HD2";
		atoms[15].pdb_type = "HE1";
		atoms[16].pdb_type = "HE2";
		atoms[17].pdb_type = "HZ";
		atoms[18].pdb_type = "C";
		atoms[19].pdb_type = "O";
	}
	//Atom aromatic_plane;
};
static factory::Registrar<F4> F4_dummy("F");

class G2 : public Residue
{
public:

	G2(){std::vector<Atom> dummy(7,Atom()); atoms = dummy;main_rotatory_bonds = 2;side_rotatory_bonds=0;resname = "GLY";res1 = "G";}
	Residue* Spawn(){return new G2;}

	Atom* GenAll(Atom * junction_atom)
	{
		atoms[0] = N(junction_atom,180.d);	// peptide nitrogen has dihedral 180, peptide oxygen has dihedral 0 (arbitrary choice)
		atoms[1] = H(&atoms[0],180.d);	// trans to oxygen
		atoms[2] = CA(&atoms[0],0.d);	main_rotatory_atoms.push_back(&atoms[2]); // trans
		atoms[3] = H(&atoms[2],60.d);
		atoms[4] = H(&atoms[2],180.d);

		atoms[5] = C(&atoms[2],300.d);	atoms[5].plane_def = &atoms[6];	main_rotatory_atoms.push_back(&atoms[5]);
		atoms[6] = O(&atoms[5],0.d);

		junction = &(*----atoms.end()); return junction;
	}
	void NameMe(){
		atoms[0].pdb_type = "N";
		atoms[1].pdb_type = "H";
		atoms[2].pdb_type = "CA";
		atoms[3].pdb_type = "1HA";
		atoms[4].pdb_type = "2HA";
		atoms[5].pdb_type = "C";
		atoms[6].pdb_type = "O";
	}
};
static factory::Registrar<G2> G2_dummy("G");

class H4 : public Residue
{
public:

	H4(){std::vector<Atom> dummy(17,Atom()); atoms = dummy;main_rotatory_bonds = 2;side_rotatory_bonds=2;resname = "HIS";res1 = "H";}
	Residue* Spawn(){return new H4;}

	Atom* GenAll(Atom * junction_atom)
	{
		atoms[0] = N(junction_atom,180.d);	// peptide nitrogen has dihedral 180, peptide oxygen has dihedral 0 (arbitrary choice)
		atoms[1] = H(&atoms[0],180.d);	// trans to oxygen
		atoms[2] = CA(&atoms[0],0.d);	main_rotatory_atoms.push_back(&atoms[2]); // trans
		atoms[3] = H(&atoms[2],180.d);

		atoms[4] = CA(&atoms[2],300.d);	side_rotatory_atoms.push_back(&atoms[4]);	//CB
		atoms[5] = H(&atoms[4],0.d);
		atoms[6] = H(&atoms[4],120.d);

		atoms[7] = C(&atoms[4],240.d);	atoms[7].plane_def = atoms[7].parent;	atoms[7].angle = (126 * M_PI / 180.d);	side_rotatory_atoms.push_back(&atoms[7]);		//CG
		atoms[8] = N(&atoms[7],0.d);		atoms[8].plane_def = atoms[8].parent;	atoms[8].angle = (117 * M_PI / 180.d);		//ND
		atoms[9] = C(&atoms[7],180.d);	atoms[9].plane_def = atoms[9].parent;	atoms[9].angle = (117 * M_PI / 180);		//CD

		atoms[10] = C(&atoms[8],0.d);	atoms[10].angle = (117 * M_PI / 180.d);;	//CE
		atoms[11] = N(&atoms[9],0.d);	atoms[11].angle = (117 * M_PI / 180.d);;	//NE

		atoms[10].children.push_back(&atoms[11]);
		atoms[11].children.push_back(&atoms[10]);

		atoms[12] = H(&atoms[9],180.d);	// hydrogen
		atoms[13] = H(&atoms[10],0.d);	// hydrogen
		atoms[14] = H(&atoms[8],0.d);	// delta hydrogen
		atoms[14] = H(&atoms[11],0.d);	// epsilon hydrogen
		
		atoms[15] = C(&atoms[2],60.d);	atoms[15].plane_def = &atoms[16];	main_rotatory_atoms.push_back(&atoms[15]);
		atoms[16] = O(&atoms[15],0.d);

		junction = &(*----atoms.end()); return junction;
	}
	void NameMe(){
		atoms[0].pdb_type = "N";
		atoms[1].pdb_type = "H";
		atoms[2].pdb_type = "CA";
		atoms[3].pdb_type = "HA";
		atoms[4].pdb_type = "CB";
		atoms[5].pdb_type = "1HB";
		atoms[6].pdb_type = "2HB";
		atoms[7].pdb_type = "CG";
		atoms[8].pdb_type = "ND1";
		atoms[9].pdb_type = "CD2";
		atoms[10].pdb_type = "CE1";
		atoms[11].pdb_type = "NE2";
		atoms[12].pdb_type = "HD2";
		atoms[13].pdb_type = "HE1";
		atoms[14].pdb_type = "HE2";
		atoms[15].pdb_type = "C";
		atoms[16].pdb_type = "O";
	}
};
static factory::Registrar<H4> H4_dummy("H");

class I4 : public Residue
{
public:

	I4(){std::vector<Atom> dummy(19,Atom()); atoms = dummy;main_rotatory_bonds = 2;side_rotatory_bonds=2;resname = "ILE";res1 = "I";}
	Residue* Spawn(){return new I4;}

	Atom* GenAll(Atom * junction_atom)
	{
		atoms[0] = N(junction_atom,180.d);	// peptide nitrogen has dihedral 180, peptide oxygen has dihedral 0 (arbitrary choice)
		atoms[1] = H(&atoms[0],180.d);	// trans to oxygen
		atoms[2] = CA(&atoms[0],0.d);	main_rotatory_atoms.push_back(&atoms[2]); // trans
		atoms[3] = H(&atoms[2],180.d);

		atoms[4] = CA(&atoms[2],300.d);	side_rotatory_atoms.push_back(&atoms[4]);	//CB
		atoms[5] = H(&atoms[4],0.d);

		atoms[6] = CA(&atoms[4],240.d);	side_rotatory_atoms.push_back(&atoms[6]);	//CG
		atoms[7] = H(&atoms[6],0.d);
		atoms[8] = H(&atoms[6],120.d);

		atoms[9] = CA(&atoms[6],240.d);	//CD
		atoms[10] = H(&atoms[9],0.d);
		atoms[11] = H(&atoms[9],120.d);
		atoms[12] = H(&atoms[9],240.d);

		atoms[13] = CA(&atoms[4],120.d); //side-methyl
		atoms[14] = H(&atoms[13],0.d);
		atoms[15] = H(&atoms[13],120.d);
		atoms[16] = H(&atoms[13],240.d);

		atoms[17] = C(&atoms[2],60.d);	atoms[17].plane_def = &atoms[18];	main_rotatory_atoms.push_back(&atoms[17]);
		atoms[18] = O(&atoms[17],0.d);

		junction = &(*----atoms.end()); return junction;
	}
	void NameMe(){
		atoms[0].pdb_type = "N";
		atoms[1].pdb_type = "H";
		atoms[2].pdb_type = "CA";
		atoms[3].pdb_type = "HA";
		atoms[4].pdb_type = "CB";
		atoms[5].pdb_type = "1HB";
		atoms[6].pdb_type = "CG1";
		atoms[7].pdb_type = "1HG1";
		atoms[8].pdb_type = "2HG1";
		atoms[9].pdb_type = "CD1";
		atoms[10].pdb_type = "1HD1";
		atoms[11].pdb_type = "2HD1";
		atoms[12].pdb_type = "3HD1";
		atoms[13].pdb_type = "CG2";
		atoms[14].pdb_type = "1HG2";
		atoms[15].pdb_type = "2HG2";
		atoms[16].pdb_type = "3HG2";
		atoms[17].pdb_type = "C";
		atoms[18].pdb_type = "O";
	}
};
static factory::Registrar<I4> I4_dummy("I");

class K6 : public Residue
{
public:

	K6(){std::vector<Atom> dummy(22,Atom()); atoms = dummy;main_rotatory_bonds = 2;side_rotatory_bonds=4;resname = "LYS";res1 = "K";}
	Residue* Spawn(){return new K6;}

	Atom* GenAll(Atom * junction_atom)
	{
		atoms[0] = N(junction_atom,180.d);	// peptide nitrogen has dihedral 180, peptide oxygen has dihedral 0 (arbitrary choice)
		atoms[1] = H(&atoms[0],180.d);	// trans to oxygen
		atoms[2] = CA(&atoms[0],0.d);	main_rotatory_atoms.push_back(&atoms[2]); // trans
		atoms[3] = H(&atoms[2],180.d);

		atoms[4] = CA(&atoms[2],300.d);	side_rotatory_atoms.push_back(&atoms[4]);	//CB
		//atoms[19] = C(&atoms[2],300.d);
		//atoms[20] = O(&atoms[19],0
		atoms[5] = H(&atoms[4],0.d);
		atoms[6] = H(&atoms[4],120.d);

		atoms[7] = CA(&atoms[4],240.d);	side_rotatory_atoms.push_back(&atoms[7]);	//CG
		atoms[8] = H(&atoms[7],0.d);
		atoms[9] = H(&atoms[7],120.d);

		atoms[10] = CA(&atoms[7],240.d);	side_rotatory_atoms.push_back(&atoms[10]);	//CD
		atoms[11] = H(&atoms[10],0.d);
		atoms[12] = H(&atoms[10],120.d);

		atoms[13] = CA(&atoms[10],240.d);	side_rotatory_atoms.push_back(&atoms[13]);	//CE
		atoms[14] = H(&atoms[13],0.d);
		atoms[15] = H(&atoms[13],120.d);

		atoms[16] = NZ(&atoms[13],240.d);	//NZ
		atoms[17] = H(&atoms[16],0.d);
		atoms[18] = H(&atoms[16],120.d);
		atoms[19] = H(&atoms[16],240.d);

		atoms[20] = C(&atoms[2],60.d);	atoms[20].plane_def = &atoms[21];	main_rotatory_atoms.push_back(&atoms[20]);
		atoms[21] = O(&atoms[20],0.d);

		junction = &(*----atoms.end()); return junction;
	}
	void NameMe(){
		atoms[0].pdb_type = "N";
		atoms[1].pdb_type = "H";
		atoms[2].pdb_type = "CA";
		atoms[3].pdb_type = "HA";
		atoms[4].pdb_type = "CB";
		atoms[5].pdb_type = "1HB";
		atoms[6].pdb_type = "2HB";
		atoms[7].pdb_type = "CG";
		atoms[8].pdb_type = "1HG";
		atoms[9].pdb_type = "1HG";
		atoms[10].pdb_type = "CD";
		atoms[11].pdb_type = "1HD";
		atoms[12].pdb_type = "2HD";
		atoms[13].pdb_type = "CE";
		atoms[14].pdb_type = "1HE";
		atoms[15].pdb_type = "2HE";
		atoms[16].pdb_type = "NZ";
		atoms[17].pdb_type = "1HZ";
		atoms[18].pdb_type = "2HZ";
		atoms[19].pdb_type = "3HZ";
		atoms[20].pdb_type = "C";
		atoms[21].pdb_type = "O";
	}
};
static factory::Registrar<K6> K6_dummy("K");

class L4 : public Residue
{
public:

	L4(){std::vector<Atom> dummy(19,Atom()); atoms = dummy;main_rotatory_bonds = 2;side_rotatory_bonds=2;resname = "LEU";res1 = "L";}
	Residue* Spawn(){return new L4;}

	Atom* GenAll(Atom * junction_atom)
	{
		atoms[0] = N(junction_atom,180.d);	// peptide nitrogen has dihedral 180, peptide oxygen has dihedral 0 (arbitrary choice)
		atoms[1] = H(&atoms[0],180.d);	// trans to oxygen
		atoms[2] = CA(&atoms[0],0.d);	main_rotatory_atoms.push_back(&atoms[2]); // trans
		atoms[3] = H(&atoms[2],180.d);

		atoms[4] = CA(&atoms[2],300.d);	side_rotatory_atoms.push_back(&atoms[4]);	//CB
		atoms[5] = H(&atoms[4],0.d);
		atoms[6] = H(&atoms[4],120.d);

		atoms[7] = CA(&atoms[4],240.d);	side_rotatory_atoms.push_back(&atoms[7]);	//CG
		atoms[8] = H(&atoms[7],0.d);

		atoms[9] = CA(&atoms[7],120.d);	//CD
		atoms[10] = H(&atoms[9],0.d);
		atoms[11] = H(&atoms[9],120.d);
		atoms[12] = H(&atoms[9],240.d);

		atoms[13] = CA(&atoms[7],240.d); //side-methyl
		atoms[14] = H(&atoms[13],0.d);
		atoms[15] = H(&atoms[13],120.d);
		atoms[16] = H(&atoms[13],240.d);

		atoms[17] = C(&atoms[2],60.d);	atoms[17].plane_def = &atoms[18];	main_rotatory_atoms.push_back(&atoms[17]);
		atoms[18] = O(&atoms[17],0.d);

		junction = &(*----atoms.end()); return junction;
	}
	void NameMe(){
		atoms[0].pdb_type = "N";
		atoms[1].pdb_type = "H";
		atoms[2].pdb_type = "CA";
		atoms[3].pdb_type = "HA";
		atoms[4].pdb_type = "CB";
		atoms[5].pdb_type = "1HB";
		atoms[6].pdb_type = "2HB";
		atoms[7].pdb_type = "CG";
		atoms[8].pdb_type = "1HG";
		atoms[9].pdb_type = "CD1";
		atoms[10].pdb_type = "1HD1";
		atoms[11].pdb_type = "2HD1";
		atoms[12].pdb_type = "3HD1";
		atoms[13].pdb_type = "CD2";
		atoms[14].pdb_type = "1HD2";
		atoms[15].pdb_type = "2HD2";
		atoms[16].pdb_type = "3HD2";
		atoms[17].pdb_type = "C";
		atoms[18].pdb_type = "O";
	}
};
static factory::Registrar<L4> L4_dummy("L");

class M5 : public Residue
{
public:

	M5(){std::vector<Atom> dummy(17,Atom()); atoms = dummy;main_rotatory_bonds = 2;side_rotatory_bonds=3;resname = "MET";res1 = "M";}
	Residue* Spawn(){return new M5;}

	Atom* GenAll(Atom * junction_atom)
	{
		atoms[0] = N(junction_atom,180.d);	// peptide nitrogen has dihedral 180, peptide oxygen has dihedral 0 (arbitrary choice)
		atoms[1] = H(&atoms[0],180.d);	// trans to oxygen
		atoms[2] = CA(&atoms[0],0.d);	main_rotatory_atoms.push_back(&atoms[2]); // trans
		atoms[3] = H(&atoms[2],180.d);

		atoms[4] = CA(&atoms[2],300.d);	side_rotatory_atoms.push_back(&atoms[4]);	//CB
		atoms[5] = H(&atoms[4],0.d);
		atoms[6] = H(&atoms[4],120.d);

		atoms[7] = CA(&atoms[4],240.d);	side_rotatory_atoms.push_back(&atoms[7]);	//CG
		atoms[8] = H(&atoms[7],0.d);
		atoms[9] = H(&atoms[7],120.d);

		atoms[10] = SM(&atoms[7],240.d);	side_rotatory_atoms.push_back(&atoms[10]);	//SD

		atoms[11] = CA(&atoms[10],0.d);		//CE
		atoms[12] = H(&atoms[11],0.d);
		atoms[13] = H(&atoms[11],120.d);
		atoms[14] = H(&atoms[11],240.d);

		atoms[15] = C(&atoms[2],60.d);	atoms[15].plane_def = &atoms[16];	main_rotatory_atoms.push_back(&atoms[15]);
		atoms[16] = O(&atoms[15],0.d);

		junction = &(*----atoms.end()); return junction;
	}
	void NameMe(){
		atoms[0].pdb_type = "N";
		atoms[1].pdb_type = "H";
		atoms[2].pdb_type = "CA";
		atoms[3].pdb_type = "HA";
		atoms[4].pdb_type = "CB";
		atoms[5].pdb_type = "1HB";
		atoms[6].pdb_type = "2HB";
		atoms[7].pdb_type = "CG";
		atoms[8].pdb_type = "1HG";
		atoms[9].pdb_type = "2HG";
		atoms[10].pdb_type = "SD";
		atoms[11].pdb_type = "CE";
		atoms[12].pdb_type = "1HE";
		atoms[13].pdb_type = "2HE";
		atoms[14].pdb_type = "3HE";
		atoms[15].pdb_type = "C";
		atoms[16].pdb_type = "O";
	}
};
static factory::Registrar<M5> M5_dummy("M");

class N4 : public Residue
{
public:

	N4(){std::vector<Atom> dummy(14,Atom()); atoms = dummy;main_rotatory_bonds = 2;side_rotatory_bonds = 2;resname = "ASN";res1 = "N";}
	Residue* Spawn(){return new N4;}

	Atom* GenAll(Atom * junction_atom)
	{
		atoms[0] = N(junction_atom,180.d);	// peptide nitrogen has dihedral 180, peptide oxygen has dihedral 0 (arbitrary choice)
		atoms[1] = H(&atoms[0],180.d);	// trans to oxygen
		atoms[2] = CA(&atoms[0],0.d);	main_rotatory_atoms.push_back(&atoms[2]); // trans
		atoms[3] = H(&atoms[2],180.d);

		atoms[4] = CA(&atoms[2],300.d);	side_rotatory_atoms.push_back(&atoms[4]);	//CB
		atoms[5] = H(&atoms[4],0.d);
		atoms[6] = H(&atoms[4],120.d);

		atoms[7] = C(&atoms[4],240.d);	side_rotatory_atoms.push_back(&atoms[7]);	//CD
		atoms[8] = O(&atoms[7],0.d); atoms[7].plane_def = &atoms[8];
		atoms[9] = N(&atoms[7],180.d);
		atoms[10] = H(&atoms[9],0.d);
		atoms[11] = H(&atoms[9],180.d);

		atoms[12] = C(&atoms[2],60.d);	atoms[12].plane_def = &atoms[13];	main_rotatory_atoms.push_back(&atoms[12]);
		atoms[13] = O(&atoms[12],0.d);

		junction = &(*----atoms.end()); return junction;
	}
	void NameMe(){
		atoms[0].pdb_type = "N";
		atoms[1].pdb_type = "H";
		atoms[2].pdb_type = "CA";
		atoms[3].pdb_type = "HA";
		atoms[4].pdb_type = "CB";
		atoms[5].pdb_type = "1HB";
		atoms[6].pdb_type = "2HB";
		atoms[7].pdb_type = "CG";
		atoms[8].pdb_type = "OD";
		atoms[9].pdb_type = "ND";
		atoms[10].pdb_type = "1HN";
		atoms[11].pdb_type = "2HN";
		atoms[12].pdb_type = "C";
		atoms[13].pdb_type = "O";
	}
};
static factory::Registrar<N4> N4_dummy("N");

class O5 : public Residue
{
public:

	O5(){std::vector<Atom> dummy(19,Atom()); atoms = dummy;main_rotatory_bonds = 2;side_rotatory_bonds = 3;resname = "ORN";res1 = "O";}
	Residue* Spawn(){return new O5;}

	Atom* GenAll(Atom * junction_atom)
	{
		atoms[0] = N(junction_atom,180.d);	// peptide nitrogen has dihedral 180, peptide oxygen has dihedral 0 (arbitrary choice)
		atoms[1] = H(&atoms[0],180.d);	// trans to oxygen
		atoms[2] = CA(&atoms[0],0.d);	main_rotatory_atoms.push_back(&atoms[2]); // trans
		atoms[3] = H(&atoms[2],180.d);

		atoms[4] = CA(&atoms[2],300.d);	side_rotatory_atoms.push_back(&atoms[4]);	//CB
		//atoms[19] = C(&atoms[2],300.d);
		//atoms[20] = O(&atoms[19],0
		atoms[5] = H(&atoms[4],0.d);
		atoms[6] = H(&atoms[4],120.d);

		atoms[7] = CA(&atoms[4],240.d);	side_rotatory_atoms.push_back(&atoms[7]);	//CG
		atoms[8] = H(&atoms[7],0.d);
		atoms[9] = H(&atoms[7],120.d);

		atoms[10] = CA(&atoms[7],240.d);	side_rotatory_atoms.push_back(&atoms[10]);	//CD
		atoms[11] = H(&atoms[10],0.d);
		atoms[12] = H(&atoms[10],120.d);

		atoms[13] = NZ(&atoms[10],240.d);	//NZ
		atoms[14] = H(&atoms[13],0.d);
		atoms[15] = H(&atoms[13],120.d);
		atoms[16] = H(&atoms[13],240.d);

		atoms[17] = C(&atoms[2],60.d);	atoms[17].plane_def = &atoms[18];	main_rotatory_atoms.push_back(&atoms[17]);
		atoms[18] = O(&atoms[17],0.d);

		junction = &(*----atoms.end()); return junction;
	}
	void NameMe(){
		atoms[0].pdb_type = "N";
		atoms[1].pdb_type = "H";
		atoms[2].pdb_type = "CA";
		atoms[3].pdb_type = "HA";
		atoms[4].pdb_type = "CB";
		atoms[5].pdb_type = "1HB";
		atoms[6].pdb_type = "2HB";
		atoms[7].pdb_type = "CG";
		atoms[8].pdb_type = "1HG";
		atoms[9].pdb_type = "1HG";
		atoms[10].pdb_type = "CD";
		atoms[11].pdb_type = "1HD";
		atoms[12].pdb_type = "2HD";
		atoms[13].pdb_type = "NE";
		atoms[14].pdb_type = "1HE";
		atoms[15].pdb_type = "2HE";
		atoms[16].pdb_type = "3HE";
		atoms[17].pdb_type = "C";
		atoms[18].pdb_type = "O";
	}
};
static factory::Registrar<O5> O5_dummy("O");

class P1 : public Residue
{
public:

	P1(){std::vector<Atom> dummy(14,Atom()); atoms = dummy;main_rotatory_bonds = 1;side_rotatory_bonds=0;resname = "PRO";res1 = "P";}
	Residue* Spawn(){return new P1;}

	Atom* GenAll(Atom * junction_atom)
	{
		atoms[0] = N(junction_atom,180.d);	// peptide nitrogen has dihedral 180, peptide oxygen has dihedral 0 (arbitrary choice)
		atoms[1] = CA(&atoms[0],0.d);	// trans
		atoms[2] = H(&atoms[1],220.d);

		atoms[3] = CA(&atoms[1],20.d);	atoms[3].angle = (105 * M_PI / 180.d); atoms[3].bond_length = 1.45; //CB
		atoms[4] = H(&atoms[3],300.d);
		atoms[5] = H(&atoms[3],60.d);

		atoms[6] = CA(&atoms[3],160.d);	atoms[6].angle = (105 * M_PI / 180.d); atoms[6].bond_length = 1.45; //CG
		atoms[7] = H(&atoms[6],60.d);
		atoms[8] = H(&atoms[6],300.d);

		atoms[9] = CA(&atoms[6],200.d); atoms[9].angle = (105 * M_PI / 180.d); atoms[9].bond_length = 1.45;	//CD // trans to oxygen
		atoms[10] = H(&atoms[9],60.d);
		atoms[11] = H(&atoms[9],300.d);

//		atoms[9].children.push_back(&atoms[0]);
//		atoms[0].children.push_back(&atoms[9]);

		atoms[12] = C(&atoms[1],120.d);	atoms[12].plane_def = &atoms[13];	main_rotatory_atoms.push_back(&atoms[12]);
		atoms[13] = O(&atoms[12],0.d);

		junction = &(*----atoms.end()); return junction;
	}
	void NameMe(){
		atoms[0].pdb_type = "N";
		atoms[1].pdb_type = "CA";
		atoms[2].pdb_type = "HA";
		atoms[3].pdb_type = "CB";
		atoms[4].pdb_type = "1HB";
		atoms[5].pdb_type = "2HB";
		atoms[6].pdb_type = "CG";
		atoms[7].pdb_type = "1HG";
		atoms[8].pdb_type = "1HG";
		atoms[9].pdb_type = "CD";
		atoms[10].pdb_type = "1HD";
		atoms[11].pdb_type = "2HD";
		atoms[12].pdb_type = "C";
		atoms[13].pdb_type = "O";
	}
};
static factory::Registrar<P1> P1_dummy("P");

class Q5 : public Residue
{
public:

	Q5(){std::vector<Atom> dummy(17,Atom()); atoms = dummy;main_rotatory_bonds = 2;side_rotatory_bonds = 3;resname = "GLN";res1 = "Q";}
	Residue* Spawn(){return new Q5;}

	Atom* GenAll(Atom * junction_atom)
	{
		atoms[0] = N(junction_atom,180.d);	// peptide nitrogen has dihedral 180, peptide oxygen has dihedral 0 (arbitrary choice)
		atoms[1] = H(&atoms[0],180.d);	// trans to oxygen
		atoms[2] = CA(&atoms[0],0.d);	main_rotatory_atoms.push_back(&atoms[2]); // trans
		atoms[3] = H(&atoms[2],180.d);

		atoms[4] = CA(&atoms[2],300.d);	side_rotatory_atoms.push_back(&atoms[4]);	//CB
		atoms[5] = H(&atoms[4],0.d);
		atoms[6] = H(&atoms[4],120.d);

		atoms[7] = CA(&atoms[4],240.d);	side_rotatory_atoms.push_back(&atoms[7]);	//CG
		atoms[8] = H(&atoms[7],0.d);
		atoms[9] = H(&atoms[7],120.d);

		atoms[10] = C(&atoms[7],240.d);	side_rotatory_atoms.push_back(&atoms[10]);	//CD
		atoms[11] = O(&atoms[10],0.d); atoms[10].plane_def = &atoms[11];
		atoms[12] = N(&atoms[10],180.d);
		atoms[13] = H(&atoms[12],0.d);
		atoms[14] = H(&atoms[12],180.d);

		atoms[15] = C(&atoms[2],60.d);	atoms[15].plane_def = &atoms[16];	main_rotatory_atoms.push_back(&atoms[15]);
		atoms[16] = O(&atoms[15],0.d);

		junction = &(*----atoms.end()); return junction;
	}
	void NameMe(){
		atoms[0].pdb_type = "N";
		atoms[1].pdb_type = "H";
		atoms[2].pdb_type = "CA";
		atoms[3].pdb_type = "HA";
		atoms[4].pdb_type = "CB";
		atoms[5].pdb_type = "1HB";
		atoms[6].pdb_type = "2HB";
		atoms[7].pdb_type = "CG";
		atoms[8].pdb_type = "1HG";
		atoms[9].pdb_type = "2HG";
		atoms[10].pdb_type = "CD";
		atoms[11].pdb_type = "OE";
		atoms[12].pdb_type = "NE";
		atoms[13].pdb_type = "1HN";
		atoms[14].pdb_type = "2HN";
		atoms[15].pdb_type = "C";
		atoms[16].pdb_type = "O";
	}
};
static factory::Registrar<Q5> Q5_dummy("Q");

class R6: public Residue
{
public:

	R6(){std::vector<Atom> dummy(24,Atom()); atoms = dummy;main_rotatory_bonds = 2;side_rotatory_bonds = 4;resname = "ARG";res1 = "R";}
	Residue* Spawn(){return new R6;}

	Atom* GenAll(Atom * junction_atom)
	{
		atoms[0] = N(junction_atom,180.d);	// peptide nitrogen has dihedral 180, peptide oxygen has dihedral 0 (arbitrary choice)
		atoms[1] = H(&atoms[0],180.d);	// trans to oxygen
		atoms[2] = CA(&atoms[0],0.d);	main_rotatory_atoms.push_back(&atoms[2]); // trans
		atoms[3] = H(&atoms[2],180.d);

		atoms[4] = CA(&atoms[2],300.d);	side_rotatory_atoms.push_back(&atoms[4]);	//CB
		atoms[5] = H(&atoms[4],0.d);
		atoms[6] = H(&atoms[4],120.d);

		atoms[7] = CA(&atoms[4],240.d);	side_rotatory_atoms.push_back(&atoms[7]);	//CG
		atoms[8] = H(&atoms[7],0.d);
		atoms[9] = H(&atoms[7],120.d);

		atoms[10] = CA(&atoms[7],240.d);	side_rotatory_atoms.push_back(&atoms[10]);	//CD
		atoms[11] = H(&atoms[10],0.d);
		atoms[12] = H(&atoms[10],120.d);

		atoms[13] = N(&atoms[10],240.d); atoms[13].plane_def = &atoms[14]; side_rotatory_atoms.push_back(&atoms[13]);
		atoms[14] = H(&atoms[13],180.d);

		atoms[15] = C(&atoms[13],0.d);	//CZ

		atoms[16] = N(&atoms[15],0.d);
		atoms[17] = H(&atoms[16],0.d);
		atoms[18] = H(&atoms[16],180.d);

		atoms[19] = N(&atoms[15],180.d);
		atoms[20] = H(&atoms[19],0.d);
		atoms[21] = H(&atoms[19],180.d);

		atoms[22] = C(&atoms[2],60.d);	atoms[22].plane_def = &atoms[23];	main_rotatory_atoms.push_back(&atoms[22]);
		atoms[23] = O(&atoms[22],0.d);

		junction = &(*----atoms.end()); return junction;
	}
	void NameMe(){
		atoms[0].pdb_type = "N";
		atoms[1].pdb_type = "H";
		atoms[2].pdb_type = "CA";
		atoms[3].pdb_type = "HA";
		atoms[4].pdb_type = "CB";
		atoms[5].pdb_type = "1HB";
		atoms[6].pdb_type = "2HB";
		atoms[7].pdb_type = "CG";
		atoms[8].pdb_type = "1HG";
		atoms[9].pdb_type = "2HG";
		atoms[10].pdb_type = "CD";
		atoms[11].pdb_type = "1HD";
		atoms[12].pdb_type = "2HD";
		atoms[13].pdb_type = "NE";
		atoms[14].pdb_type = "HN";
		atoms[15].pdb_type = "CZ";
		atoms[16].pdb_type = "NH1";
		atoms[17].pdb_type = "1HH1";
		atoms[18].pdb_type = "2HH1";
		atoms[19].pdb_type = "NH1";
		atoms[20].pdb_type = "1HH2";
		atoms[21].pdb_type = "2HH2";
		atoms[22].pdb_type = "C";
		atoms[23].pdb_type = "O";
	}
};
static factory::Registrar<R6> R6_dummy("R");

class S3 : public Residue
{
public:

	S3(){std::vector<Atom> dummy(11,Atom()); atoms = dummy;main_rotatory_bonds = 2;side_rotatory_bonds=1;resname = "SER";res1 = "S";}
	Residue* Spawn(){return new S3;}

		Atom* GenAll(Atom * junction_atom)
	{
		atoms[0] = N(junction_atom,180.d);	// peptide nitrogen has dihedral 180, peptide oxygen has dihedral 0 (arbitrary choice)
		atoms[1] = H(&atoms[0],180.d);	// trans to oxygen
		atoms[2] = CA(&atoms[0],0.d);	main_rotatory_atoms.push_back(&atoms[2]); // trans
		atoms[3] = H(&atoms[2],180.d);

		atoms[4] = CA(&atoms[2],300.d);	side_rotatory_atoms.push_back(&atoms[4]); //CB
		atoms[5] = H(&atoms[4],0.d);
		atoms[6] = H(&atoms[4],120.d);
		
		atoms[7] = O(&atoms[4],240.d);	//OG
		atoms[8] = H(&atoms[7],0.d);

		atoms[9] = C(&atoms[2],60.d);	atoms[9].plane_def = &atoms[10];	main_rotatory_atoms.push_back(&atoms[9]); 
		atoms[10] = O(&atoms[9],0.d);

		junction = &(*----atoms.end()); return junction;
	}
	void NameMe(){
		atoms[0].pdb_type = "N";
		atoms[1].pdb_type = "H";
		atoms[2].pdb_type = "CA";
		atoms[3].pdb_type = "HA";
		atoms[4].pdb_type = "CB";
		atoms[5].pdb_type = "1HB";
		atoms[6].pdb_type = "2HB";
		atoms[7].pdb_type = "OG";
		atoms[8].pdb_type = "HG";
		atoms[9].pdb_type = "C";
		atoms[10].pdb_type = "O";
	}
};
static factory::Registrar<S3> S3_dummy("S");

class T3 : public Residue
{
public:

	T3(){std::vector<Atom> dummy(14,Atom()); atoms = dummy;main_rotatory_bonds = 2;side_rotatory_bonds=1;resname = "THR";res1 = "T";}
	Residue* Spawn(){return new T3;}

		Atom* GenAll(Atom * junction_atom)
	{
		atoms[0] = N(junction_atom,180.d);	// peptide nitrogen has dihedral 180, peptide oxygen has dihedral 0 (arbitrary choice)
		atoms[1] = H(&atoms[0],180.d);	// trans to oxygen
		atoms[2] = CA(&atoms[0],0.d);	main_rotatory_atoms.push_back(&atoms[2]); // trans
		atoms[3] = H(&atoms[2],180.d);

		atoms[4] = CA(&atoms[2],300.d);	side_rotatory_atoms.push_back(&atoms[4]); //CB
		atoms[5] = H(&atoms[4],0.d);
		
		atoms[6] = O(&atoms[4],240.d);	//OG
		atoms[7] = H(&atoms[6],0.d);

		atoms[8] = CA(&atoms[4],120.d);
		atoms[9] = H(&atoms[8],0.d);
		atoms[10] = H(&atoms[8],120.d);
		atoms[11] = H(&atoms[8],240.d);

		atoms[12] = C(&atoms[2],60.d);	atoms[12].plane_def = &atoms[13];	main_rotatory_atoms.push_back(&atoms[12]); 
		atoms[13] = O(&atoms[12],0.d);

		junction = &(*----atoms.end()); return junction;
	}
	void NameMe(){
		atoms[0].pdb_type = "N";
		atoms[1].pdb_type = "H";
		atoms[2].pdb_type = "CA";
		atoms[3].pdb_type = "HA";
		atoms[4].pdb_type = "CB";
		atoms[5].pdb_type = "1HB";
		atoms[6].pdb_type = "OG1";
		atoms[7].pdb_type = "HO";
		atoms[8].pdb_type = "CG1";
		atoms[9].pdb_type = "1HG";
		atoms[10].pdb_type = "2HG";
		atoms[11].pdb_type = "3HG";
		atoms[12].pdb_type = "C";
		atoms[13].pdb_type = "O";
	}
};
static factory::Registrar<T3> T3_dummy("T");

class V3 : public Residue
{
public:

	V3(){std::vector<Atom> dummy(16,Atom()); atoms = dummy;main_rotatory_bonds = 2;side_rotatory_bonds = 1;resname = "VAL";res1 = "V";}
	Residue* Spawn(){return new V3;}

	Atom* GenAll(Atom * junction_atom)
	{
		atoms[0] = N(junction_atom,180.d);	// peptide nitrogen has dihedral 180, peptide oxygen has dihedral 0 (arbitrary choice)
		atoms[1] = H(&atoms[0],180.d);	// trans to oxygen
		atoms[2] = CA(&atoms[0],0.d);	main_rotatory_atoms.push_back(&atoms[2]); // trans
		atoms[3] = H(&atoms[2],180.d);

		atoms[4] = CA(&atoms[2],300.d);	side_rotatory_atoms.push_back(&atoms[4]);	//CB
		atoms[5] = H(&atoms[4],0.d);

		atoms[6] = CA(&atoms[4],120.d);	//CG1
		atoms[7] = H(&atoms[6],0.d);
		atoms[8] = H(&atoms[6],120.d);
		atoms[9] = H(&atoms[6],240.d);

		atoms[10] = CA(&atoms[4],240.d); //CG2
		atoms[11] = H(&atoms[10],0.d);
		atoms[12] = H(&atoms[10],120.d);
		atoms[13] = H(&atoms[10],240.d);

		atoms[14] = C(&atoms[2],60.d);	atoms[14].plane_def = &atoms[15];	main_rotatory_atoms.push_back(&atoms[14]);
		atoms[15] = O(&atoms[14],0.d);

		junction = &(*----atoms.end()); return junction;
	}
	void NameMe(){
		atoms[0].pdb_type = "N";
		atoms[1].pdb_type = "H";
		atoms[2].pdb_type = "CA";
		atoms[3].pdb_type = "HA";
		atoms[4].pdb_type = "CB";
		atoms[5].pdb_type = "1HB";
		atoms[6].pdb_type = "CG1";
		atoms[7].pdb_type = "1HG1";
		atoms[8].pdb_type = "2HG1";
		atoms[9].pdb_type = "3HG1";
		atoms[10].pdb_type = "CG2";
		atoms[11].pdb_type = "1HG2";
		atoms[12].pdb_type = "2HG2";
		atoms[13].pdb_type = "3HG2";
		atoms[14].pdb_type = "C";
		atoms[15].pdb_type = "O";
	}
};
static factory::Registrar<V3> V3_dummy("V");

class Y4 : public Residue
{
public:

	Y4(){std::vector<Atom> dummy(21,Atom()); atoms = dummy;main_rotatory_bonds = 2;side_rotatory_bonds = 2;resname = "TYR";res1 = "Y";}
	Residue* Spawn(){return new Y4;}

	Atom* GenAll(Atom * junction_atom)
	{
		atoms[0] = N(junction_atom,180.d);	// peptide nitrogen has dihedral 180, peptide oxygen has dihedral 0 (arbitrary choice)
		atoms[1] = H(&atoms[0],180.d);	// trans to oxygen
		atoms[2] = CA(&atoms[0],0.d);	main_rotatory_atoms.push_back(&atoms[2]); // trans
		atoms[3] = H(&atoms[2],180.d);

		atoms[4] = CA(&atoms[2],300.d);	side_rotatory_atoms.push_back(&atoms[4]);	//CB
		atoms[5] = H(&atoms[4],0.d);
		atoms[6] = H(&atoms[4],120.d);

		//aromatic_plane.coordinates = atoms[5]-atoms[4] + atoms[6]-atoms[4];

		atoms[7] = C(&atoms[4],240.d);	atoms[7].plane_def = atoms[7].parent;		side_rotatory_atoms.push_back(&atoms[7]);	//CG
		atoms[8] = C(&atoms[7],0.d);		atoms[8].plane_def = atoms[8].parent;		//CD
		atoms[9] = C(&atoms[7],180.d);	atoms[9].plane_def = atoms[9].parent;		//CD

		atoms[10] = C(&atoms[8],0.d);											//CE
		atoms[11] = C(&atoms[9],0.d);		atoms[11].plane_def = atoms[11].parent;	//CE

		atoms[12] = C(&atoms[11],180.d);										//CF

		atoms[12].children.push_back(&atoms[10]);
		atoms[10].children.push_back(&atoms[12]);

		atoms[13] = H(&atoms[8],180.d - atoms[10].dihedral);	//aromatic hydrogens
		atoms[14] = H(&atoms[9],180.d - atoms[11].dihedral);	//aromatic hydrogens
		atoms[15] = H(&atoms[10],0.d);	//aromatic hydrogens
		atoms[16] = H(&atoms[11],180.d - atoms[12].dihedral);	//aromatic hydrogens
		atoms[17] = O(&atoms[12],0.d);	//aromatic hydrogens
		atoms[18] = H(&atoms[17],0.d);	//aromatic hydrogens

		
		
		atoms[19] = C(&atoms[2],60.d);	atoms[19].plane_def = &atoms[20];	main_rotatory_atoms.push_back(&atoms[19]);
		atoms[20] = O(&atoms[19],0.d);

		junction = &(*----atoms.end()); return junction;
	}
	void NameMe(){
		atoms[0].pdb_type = "N";
		atoms[1].pdb_type = "H";
		atoms[2].pdb_type = "CA";
		atoms[3].pdb_type = "HA";
		atoms[4].pdb_type = "CB";
		atoms[5].pdb_type = "1HB";
		atoms[6].pdb_type = "2HB";
		atoms[7].pdb_type = "CG";
		atoms[8].pdb_type = "CD1";
		atoms[9].pdb_type = "CD2";
		atoms[10].pdb_type = "CE1";
		atoms[11].pdb_type = "CE2";
		atoms[12].pdb_type = "CZ";
		atoms[13].pdb_type = "HD1";
		atoms[14].pdb_type = "HD2";
		atoms[15].pdb_type = "HE1";
		atoms[16].pdb_type = "HE2";
		atoms[17].pdb_type = "OH";
		atoms[18].pdb_type = "HO";
		atoms[19].pdb_type = "C";
		atoms[20].pdb_type = "O";
	}
	//Atom aromatic_plane;
};
static factory::Registrar<Y4> Y4_dummy("Y");

class W4 : public Residue
{
public:

	W4(){std::vector<Atom> dummy(24,Atom()); atoms = dummy;main_rotatory_bonds = 2; side_rotatory_bonds = 2;resname = "TRP";res1 = "W";}
	Residue* Spawn(){return new W4;}

	Atom* GenAll(Atom * junction_atom)
	{
		atoms[0] = N(junction_atom,180.d);	// peptide nitrogen has dihedral 180, peptide oxygen has dihedral 0 (arbitrary choice)
		atoms[1] = H(&atoms[0],180.d);	// trans to oxygen
		atoms[2] = CA(&atoms[0],0.d);	main_rotatory_atoms.push_back(&atoms[2]); // trans
		atoms[3] = H(&atoms[2],180.d);

		atoms[4] = CA(&atoms[2],300.d);	side_rotatory_atoms.push_back(&atoms[4]);	//CB
		atoms[5] = H(&atoms[4],0.d);
		atoms[6] = H(&atoms[4],120.d);

		atoms[7] = C(&atoms[4],240.d);	atoms[7].plane_def = atoms[7].parent;	atoms[7].angle = (126 * M_PI / 180.d);	side_rotatory_atoms.push_back(&atoms[7]);		//CG
		atoms[8] = CR(&atoms[7],0.d); atoms[8].plane_def = atoms[8].parent; atoms[8].angle = (135 * M_PI / 180.d);		//CD1
		atoms[9] = CR(&atoms[7],180.d);	atoms[9].plane_def = atoms[9].parent;	atoms[9].angle = (108 * M_PI / 180);		//CD2

		atoms[10] = N(&atoms[9],0.d);	atoms[11].angle = (117 * M_PI / 180.d);;	//NE

		atoms[11] = CR(&atoms[8],180.d);	atoms[12].plane_def = atoms[12].parent;		//CG
		atoms[12] = CR(&atoms[11],0.d);	atoms[13].plane_def = atoms[13].parent;		//CG
		atoms[13] = CR(&atoms[12],180.d);	atoms[14].plane_def = atoms[14].parent;		//CG
		atoms[14] = CR(&atoms[13],180.d);	atoms[15].plane_def = atoms[15].parent;		//CG
		atoms[15] = CR(&atoms[14],180.d);	atoms[15].plane_def = atoms[15].parent;		//CG
		
		atoms[10].children.push_back(&atoms[15]);
		atoms[15].children.push_back(&atoms[10]);
		atoms[8].children.push_back(&atoms[15]);
		atoms[15].children.push_back(&atoms[8]);

		atoms[16] = H(&atoms[9],180.d);	// hydrogen
		atoms[17] = H(&atoms[12],0.d);	// hydrogen
		atoms[18] = H(&atoms[13],0.d);	// delta hydrogen
		atoms[19] = H(&atoms[14],0.d);	// epsilon hydrogen

		atoms[20] = H(&atoms[10],0.d);	// epsilon hydrogen
		atoms[21] = H(&atoms[11],180.d);	// epsilon hydrogen
		
		atoms[22] = C(&atoms[2],60.d);	atoms[22].plane_def = &atoms[23];	main_rotatory_atoms.push_back(&atoms[22]);
		atoms[23] = O(&atoms[22],0.d);

		junction = &(*----atoms.end()); return junction;
	}
	void NameMe(){
		atoms[0].pdb_type = "N";
		atoms[1].pdb_type = "H";
		atoms[2].pdb_type = "CA";
		atoms[3].pdb_type = "HA";
		atoms[4].pdb_type = "CB";
		atoms[5].pdb_type = "1HB";
		atoms[6].pdb_type = "2HB";
		atoms[7].pdb_type = "CG";
		atoms[8].pdb_type = "CD2";
		atoms[9].pdb_type = "CD1";
		atoms[10].pdb_type = "NE1";
		atoms[11].pdb_type = "CE3";
		atoms[12].pdb_type = "CZ3";
		atoms[13].pdb_type = "CH2";
		atoms[14].pdb_type = "CZ2";
		atoms[15].pdb_type = "CE2";
		atoms[16].pdb_type = "HD1";
		atoms[17].pdb_type = "HZ3";
		atoms[18].pdb_type = "HH2";
		atoms[19].pdb_type = "HZ2";
		atoms[20].pdb_type = "HE1";
		atoms[21].pdb_type = "HE2";
		atoms[22].pdb_type = "C";
		atoms[23].pdb_type = "O";
	}
};
static factory::Registrar<W4> W4_dummy("W");

class Ncap : public Residue
{
public:

	Ncap(){std::vector<Atom> dummy(6,Atom()); atoms = dummy;main_rotatory_bonds = 0;side_rotatory_bonds = 0;resname = "NAP";res1 = "";}
	Residue* Spawn(){return new Ncap;}

	Atom* GenAll(Atom * junction_atom)
	{
		atoms[0] = N(junction_atom,180.d);	// peptide nitrogen has dihedral 180, peptide oxygen has dihedral 0 (arbitrary choice)
		atoms[1] = H(&atoms[0],180.d);	// trans to oxygen
		atoms[2] = CA(&atoms[0],0.d);	// trans
		atoms[3] = H(&atoms[2],60.d);
		atoms[4] = H(&atoms[2],180.d);
		atoms[5] = H(&atoms[2],300.d);	

		junction = NULL; return junction;
	}
	void NameMe(){
		atoms[0].pdb_type = "N";
		atoms[1].pdb_type = "H";
		atoms[2].pdb_type = "CA";
		atoms[3].pdb_type = "1HA";
		atoms[4].pdb_type = "2HA";
		atoms[5].pdb_type = "3HA";
	}
};
static factory::Registrar<Ncap> Ncap_dummy("Ncap");
