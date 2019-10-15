#include <iostream>
#include "base_classes.h"
#include "factory.h"
#include <iomanip>
#include <iterator>
#include <locale>
#include <sys/stat.h>
#include <props.h>

bool file_output = false;
bool screen_output = true; 

void msg(std::string msg, std::string file = "logfile")
{
	if ( file_output )
	{
		std::ofstream myofile;
		myofile.open(file, std::fstream::out | std::fstream::app);
		myofile << msg << std::endl;
		myofile.close();
	}
	
	if ( screen_output )
	{
		std::cout << msg << std::endl;	
	}
}

long myPow(long x, long p) {
  if (p == 0) return 1;
  if (p == 1) return x;
  return x * myPow(x, p-1);
}

void subbase_to_n(long number, int base, std::string & current){
  int i = 0;
  while ( myPow(base,i+1) <= number ){i++;}
  current[i] = std::to_string( number / myPow(base,i))[0];
  if ( i != 0 ){
  	subbase_to_n(number % myPow(base,i), base, current);
  }
}

void base_to_n(long number, int base_main, int base_side, std::string & main_current, std::string & side_current)
{
	long main_number, side_number;
	side_number = (number % myPow(base_side, side_current.length()));
	main_number = (number - side_number) / myPow(base_side, side_current.length());

	subbase_to_n(main_number, base_main, main_current);
	subbase_to_n(side_number, base_side, side_current);
}

bool is_number(const std::string& s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}
//// ///////// /////

int main(){
	
	//Ccap res;
//	res.GenAll();
 msg("started");
	std::string a("Ccap");
	std::vector<Residue*> newres;
	Residue peptide;
	Atom * junction = NULL;
	double threshold_clash = 1.6;
	double hydrogen_clash = 1.6; 
	double clash_value = threshold_clash;

	newres.push_back(factory::ResidueFactory::getInstance(a));
	int res_index = 0;
	junction = newres[0]->GenAll(junction);
	newres[0]->NameMe();
	
	std::cout << "Enter value for hydrogen clash" << std::endl;
	std::cin >> a;
	hydrogen_clash = stod(a);
	std::cout << "H clash = " << hydrogen_clash << std::endl;

	std::cout << "Enter value for heavy-atom clash" << std::endl;
	std::cin >> a;
	threshold_clash = stod(a);
	std::cout << "Heavy atom clash = " << threshold_clash << std::endl;

	std::cout << "Enter number of main-chain rotations per bond" << std::endl;
	do {
		std::cin >> a;
		std::locale loc;
		if ( is_number(a)) {break;}
		std::cout << "NaN! Try again..." << std::endl;
	} while (true);

	long no_main_rotations = stoi(a);

	std::cout << "Enter number of side-chain rotations per bond" << std::endl;
	do {
		std::cin >> a;
		std::locale loc;
		if ( is_number(a)) {break;}
		std::cout << "NaN! Try again..." << std::endl;
	} while (true);

	long no_side_rotations = stoi(a);
		
	std::cout << "Enter residue starting from N-end" << std::endl;
	
	// Prepare a vector of Residues to be built
	Residue * pres;
	do {
		std::cin >> a;
		if ( a == "q" ){break;}
		pres = factory::ResidueFactory::getInstance(a);
		if ( pres == NULL){std::cout << "Unknown residue, try again..." << std::endl; continue;}
		newres.push_back(pres->Spawn());
		res_index++;
		junction = newres[res_index]->GenAll(junction);
		newres[res_index]->NameMe();
	} while (true);

	res_index++;
	newres.push_back(factory::ResidueFactory::getInstance("Ncap"));
	junction = newres[res_index]->GenAll(junction);
	newres[res_index]->NameMe();
	
	std::string dirname = "";
	for ( std::vector<Residue*>::iterator it = newres.begin(); it != newres.end(); it++ ){ dirname = dirname + (*it)->res1; };

	int errnum = mkdir(dirname.c_str(), 0744);
        if ( errnum != 0 && ! props::file_exists(dirname) ){ std::cout << "SUBMITTER: Unable to create dir " + dirname << std::endl; return 1;}

	/*/ Childrenify

	for ( std::vector<Residue*>::iterator it = newres.begin(); it != newres.end(); it++)
	{
		for ( std::vector<Atom>::iterator it2 = (*it)->atoms.begin(); it2 != (*it)->atoms.end(); it2++)
		{
			for ( std::vector<Atom>::iterator it3 = (*it)->atoms.begin(); it3 != (*it)->atoms.end(); it3++)
			{
				if ( it3->parent == &*it2 ){it2->children.push_back(&*it3);}
			}
			std::cout << (*it)->junction << std::endl;
			if ( &*it2 == (*it)->junction && it != newres.end() - 1){it2->children.push_back(&(*it)->atoms[0]);}
		}
	}
	//*/


	long total_main_rotatory_bonds = 0;
	long total_side_rotatory_bonds = 0;
	for (std::vector<Residue*>::iterator it = newres.begin(); it != newres.end(); it++){
          total_main_rotatory_bonds += (*it)->main_rotatory_bonds;
          total_side_rotatory_bonds += (*it)->side_rotatory_bonds;
        }
	
	// Prepare a vector of pointers to atoms that are rotatory
	std::vector<Atom*> final_main_rotatory_combination;
	std::vector<Atom*> final_side_rotatory_combination;

	for (std::vector<Residue*>::iterator it = newres.begin(); it != newres.end(); it++)
	{
		for(std::vector<Atom*>::iterator it_a = (*it)->main_rotatory_atoms.begin(); it_a != (*it)->main_rotatory_atoms.end(); it_a++){
			(*it_a)->rotations = no_main_rotations;
		}	
		final_main_rotatory_combination.insert(final_main_rotatory_combination.end(), (*it)->main_rotatory_atoms.begin(), (*it)->main_rotatory_atoms.end()); 
	}

	for (std::vector<Residue*>::iterator it = newres.begin(); it != newres.end(); it++)
	{
		for(std::vector<Atom*>::iterator it_a = (*it)->side_rotatory_atoms.begin(); it_a != (*it)->side_rotatory_atoms.end(); it_a++){
			(*it_a)->rotations = no_side_rotations;
		}	
		final_side_rotatory_combination.insert(final_side_rotatory_combination.end(), (*it)->side_rotatory_atoms.begin(), (*it)->side_rotatory_atoms.end()); 
	}
////////////////////


	long  total_combinations = myPow(no_main_rotations, total_main_rotatory_bonds) * myPow(no_side_rotations, total_side_rotatory_bonds); // cap rotations?
	std::cout << "Total number of conformers: " << total_combinations << std::endl;

	long not_discarded = 0;
	for (long  ii = 0; ii < total_combinations; ii++)
	{
		// Translate ii to a string of indices that define dihedrals of rotatory atoms (changing an index produces a different rotation and that's how we get conformers)
		std::string s_main (total_main_rotatory_bonds, '0'); // cap rotations?
		std::string s_side (total_side_rotatory_bonds, '0'); // cap rotations?
	
		base_to_n(ii,no_main_rotations,no_side_rotations,s_main,s_side);

//		final_rotatory_combination[0]->rotatory_index = 1; // to avoid clash of non-rotating cap
		for (int j = 0; j < total_main_rotatory_bonds - 0; j++) // start with second, end with penultimate - i.e. do not rotate caps
		{
			final_main_rotatory_combination[j]->rotatory_index = s_main[j] - '0';
		}
		for (int j = 0; j < total_side_rotatory_bonds - 0; j++) // start with second, end with penultimate - i.e. do not rotate caps
		{
			final_side_rotatory_combination[j]->rotatory_index = s_side[j] - '0';
		}

		// Generate a conformer
		int atomid = 0;
		for (int i = 0; i < newres.size(); i++)
		{
			newres[i]->resid = i;
			for (std::vector<Atom>::iterator it2 = newres[i]->atoms.begin(); it2 != newres[i]->atoms.end(); it2++)
			{
				atomid++;
				if (it2->genme){it2->GenMe();} // b_angle je nepouzivane + je false by default! always false?!
			}
		}

		// Iterate through all pairs of atoms once and look for clashes.
		//UGLY - create vector of pointers to all atoms and iterate through 2 for cycles, instead of 4

		bool clash = false;
		for (std::vector<Residue*>::iterator i = newres.begin(); i != newres.end(); i++)
		{
			for (std::vector<Atom>::iterator it = (*i)->atoms.begin(); it != (*i)->atoms.end(); it++)
			{
				if ( i == newres.end() - 1 && it == (*i)->atoms.end() - 1 ){continue;}
				std::vector<Residue*>::iterator j = i;
				std::vector<Atom>::iterator it2 = it + 1;
				if ( it + 1 == (*i)->atoms.end() ){j++; if ( j != newres.end() ){it2 = (*j)->atoms.begin();}}
				// vytvorili sme dvojicu atomov - (i,it) a (j,it2) tak, aby ten druhy bol "neskorsi" ako prvy (aby sme necheckovali 2x)	
				for (j; j != newres.end(); j++)
				{
					for (it2; it2 != (*j)->atoms.end(); it2++)
					{
						// Are it and it2 neighbours (parent-child relationship)? If so, skip proximity check
						if ( (it->parent != &(*it2)) && std::find(it->children.begin(),it->children.end(),&*it2) == it->children.end() )
						{
//							std::cout << it->type << " " << it2->type << " " << (it->coordinates - it2->coordinates).norm() << std::endl;
							// Is the distance btw them below threshold_clash? If so, declare this a clash and do not output this conformer
							if ( it->type == "H" || it2->type == "H" ){ clash_value = hydrogen_clash; } else { clash_value = threshold_clash; }
							if ( (it->coordinates - it2->coordinates).norm() < clash_value ){
								
//								std::cout << "Clash: " << std::distance(newres.begin(), i) << " " << std::distance((*i)->atoms.begin(), it) << " " << std::distance(newres.begin(), j) << " " << std::distance((*j)->atoms.begin(), it2) << " " << (it->coordinates - it2->coordinates).norm() << std::endl;
//								std::cout << "Clash: " << it->coordinates << " " << it2->coordinates << std::endl;

								clash = true; break;
							}	// there is a clash
						}
					}
					if (clash){break;}
					if ( j + 1 != newres.end() ){it2 = (*(j + 1))->atoms.begin();}
				}
				if (clash){break;}
			}
			if (clash){break;}
		}
		if (clash){std::cout << ii << " clash" << std::endl; continue;}

		std::stringstream ss;
		ss << not_discarded;

		// Print a conformer
		std::string filename = "./" + dirname + "/" + dirname + "_" + s_main + "_" + s_side;
		not_discarded++;
		std::ofstream myofile;
		myofile.open(filename);
			
		atomid = 0;
		for (int i = 0; i < newres.size(); i++)
		{
			newres[i]->resid = i;
			for (std::vector<Atom>::iterator it = newres[i]->atoms.begin(); it != newres[i]->atoms.end(); it++)
			{
				atomid++;
				it->PrintMe(myofile,atomid,newres[i]->resname,newres[i]->resid);
			}
		}

		myofile.close();

	}

	std::cout << "Done!" << std::endl;

	char c;
	std::cin >> c;
}
