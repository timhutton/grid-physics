// local:
#include "Arena.hpp"

// STL:
#include <stdexcept>
#include <algorithm>

int randInt( int n )  // [0,n)
{
    return rand() % n; // use something more uniform if needed
}

//----------------------------------------------------------------------------

const int Arena::vNx[4] = { 0, 1, 0, -1 }; // NESW
const int Arena::vNy[4] = { -1, 0, 1, 0 };

Arena::Arena(int x, int y)
    : X( x )
	, Y( y )
{
	this->grid = vector<vector<Slot>>(X,vector<Slot>(Y));
}

//----------------------------------------------------------------------------

bool Arena::isOffGrid( int x, int y ) const {
    return x < 0 || x >= X || y < 0 || y >= Y;
}

//----------------------------------------------------------------------------

bool Arena::hasAtom( int x, int y ) const {
    if( isOffGrid(x,y ) )
		throw out_of_range("Atom not on grid");

    return this->grid[x][y].has_atom;
}
    
//----------------------------------------------------------------------------

size_t Arena::addAtom( int x, int y, int type) {

    if( isOffGrid(x,y ) )
		throw out_of_range("Atom not on grid");

	Slot& slot = this->grid[x][y];
	if( slot.has_atom )
		throw exception("Grid already contains an atom at that position");

	Atom a;
	a.x = x;
	a.y = y;
    a.type = type;
	this->atoms.push_back( a );
	size_t iAtom = this->atoms.size()-1;

    slot.has_atom = true;
    slot.iAtom = iAtom;

	Group group;
	group.atoms.push_back( iAtom );
	this->groups.push_back( group );

	return iAtom;
}

//----------------------------------------------------------------------------

void Arena::makeBond( size_t a, size_t b, BondType type ) {
	if( a < 0 || a >= atoms.size() || b < 0 || b >= atoms.size() )
		throw out_of_range("Invalid atom index");
	if( a == b )
		throw invalid_argument("Cannot bond atom to itself");
    if( !isWithinFlexibleBondNeighborhood( this->atoms[a].x, this->atoms[a].y, this->atoms[b].x, this->atoms[b].y ) )
        throw invalid_argument("Atoms are too far apart to be bonded");
    if( find( begin( this->atoms[ a ].bonded_atoms ), end( this->atoms[ a ].bonded_atoms ), b ) != end( this->atoms[ a ].bonded_atoms ) )
        throw invalid_argument("Atoms are already bonded");

    switch( type ) {
        case BondType::Moore:
            this->atoms[ a ].bonded_atoms.push_back( b );
            this->atoms[ b ].bonded_atoms.push_back( a );
            break;
        case BondType::vonNeumann:
            this->atoms[ a ].rigid_bonded_atoms.push_back( b );
            this->atoms[ b ].rigid_bonded_atoms.push_back( a );
            break;
    }
}

//----------------------------------------------------------------------------

void Arena::makeGroups( FlexibilityMethod flexibility_method ) {
    switch( flexibility_method )
    {
        case FlexibilityMethod::JustAtoms: break; // no new groups to add (already have one per atom)
        case FlexibilityMethod::AllGroups: addAllGroups(); break;
        case FlexibilityMethod::Tree: makeTreeOfGroups(); break;
    }
}

//----------------------------------------------------------------------------

void Arena::makeTreeOfGroups() {
    // TODO
    /*vector<Group> working_groups = this->groups;
    while( true ) {
        vector<Group> groups_to_be_merged;
        for( const Group& group : working_groups ) {
            // if this group has exiting bonds then it needs to be merged
            if( test )
                groups_to_be_merged.push_back( group );
        }
        if( groups_to_be_merged.empty() )
            break;
        // merge the groups, 
        working_groups.clear();
        for( something ) {
            merged_group = something;
            working_groups.push_back( merged_group );
            this->groups.push_back( merged_group );
        }
    }*/
}

//----------------------------------------------------------------------------

void Arena::addAllGroups() {
    // add the rigid bonds first because they reduce the number of groups
    for( size_t iAtom = 0; iAtom < this->atoms.size(); ++iAtom ) {
        const Atom& a = this->atoms[ iAtom ];
        for( const size_t& iAtom2 : a.rigid_bonded_atoms ) {
            if( iAtom2 < iAtom ) continue; 
            addAllGroupsForNewBond( iAtom, iAtom2, BondType::vonNeumann );
        }
    }
    for( size_t iAtom = 0; iAtom < this->atoms.size(); ++iAtom ) {
        const Atom& a = this->atoms[ iAtom ];
        for( const size_t& iAtom2 : a.bonded_atoms ) {
            if( iAtom2 < iAtom ) continue; 
            addAllGroupsForNewBond( iAtom, iAtom2, BondType::Moore );
        }
    }
}

//----------------------------------------------------------------------------

void Arena::addAllGroupsForNewBond( size_t a, size_t b, BondType type ) {
	// add new groups obtained by combining pairwise every group that includes a but not b 
    // with every group that includes b but not a
	vector<Group> new_groups;
	for( const auto& ga : this->groups ) {
		if( find( begin( ga.atoms ), end( ga.atoms ), a ) == end( ga.atoms ) )
			continue;
		if( find( begin( ga.atoms ), end( ga.atoms ), b ) != end( ga.atoms ) )
			continue;
        // (ga contains a but not b)
		for( const auto& gb : this->groups ) {
			if( find( begin( gb.atoms ), end( gb.atoms ), b ) == end( gb.atoms ) )
				continue;
			if( find( begin( gb.atoms ), end( gb.atoms ), a ) != end( gb.atoms ) )
				continue;
            // (gb contains b but not a)
            // merge the two groups
			Group new_group;
            new_group.atoms.resize( ga.atoms.size() + gb.atoms.size() );
            const auto& end = set_union( ga.atoms.begin(), ga.atoms.end(), gb.atoms.begin(), gb.atoms.end(), new_group.atoms.begin() );
            new_group.atoms.resize( end - new_group.atoms.begin() );
            // add to the list if unique
            bool is_unique = true;
            for( const auto& g2 : this->groups ) {
                if( g2.atoms == new_group.atoms ) {
                    is_unique = false;
                    break;
                }
            }
			if( is_unique )
                new_groups.push_back( new_group );
        }
    }
	this->groups.insert( this->groups.end(), new_groups.begin(), new_groups.end() );
    // one optimisation: if it's a rigid bond then there's no point trying to move a and b separately
    if( type == BondType::vonNeumann ) {
        // remove any group that contains exactly one of a and b (can't move one without the other if a rigid bond)
        this->groups.erase( remove_if( begin( this->groups ), end( this->groups ), 
            GroupHasOneButNotTheOther( a, b ) ), end( this->groups ) );
    }
}

//----------------------------------------------------------------------------

bool Arena::isWithinFlexibleBondNeighborhood( int x1, int y1, int x2, int y2 ) {
    return abs( x1 - x2 ) <= 1 && abs( y1 - y2 ) <= 1;
}

//----------------------------------------------------------------------------

void Arena::update() {
    // attempt to move every group
    for( const auto& group : this->groups ) 
    {
        //const Group& group = this->groups[ rand() % this->groups.size() ];
        moveGroupRandomly( group );
    }
    doChemistry();
}

//----------------------------------------------------------------------------

void Arena::doChemistry() {
    for( int x = 0; x < this->X; ++x ) {
        for( int y = 0; y < this->Y; ++y ) {
            if( !this->grid[x][y].has_atom ) continue;
            for( int iMove = 0; iMove < 4; ++iMove ) {
                int tx = x + this->vNx[ iMove ];
                int ty = y + this->vNy[ iMove ];
                if( isOffGrid( tx, ty ) || !this->grid[tx][ty].has_atom ) continue;
                size_t iAtomA = this->grid[x][y].iAtom;
                size_t iAtomB = this->grid[tx][ty].iAtom;
                Atom& a = this->atoms[ iAtomA ];
                Atom& b = this->atoms[ iAtomB ];
                bool is_bonded = find( begin(a.bonded_atoms),end(a.bonded_atoms),iAtomB ) != end(a.bonded_atoms);
                if( !is_bonded && a.type == b.type && a.bonded_atoms.size() + b.bonded_atoms.size() < 2 ) {
                    makeBond( iAtomA, iAtomB, BondType::Moore );
                }
            }
        }
    }
}

//----------------------------------------------------------------------------

void Arena::moveGroupRandomly( const Group& group ) {
    int iMove = randInt( 4 );
    int dx = vNx[ iMove ];
    int dy = vNy[ iMove ];
    // first test: would this move stretch any bond too far?
    bool can_move = true;
    for( const size_t& iAtomIn : group.atoms ) {
        for( const size_t& iAtomOut : this->atoms[ iAtomIn ].rigid_bonded_atoms ) {
            bool b_in_group = find( begin( group.atoms ), end( group.atoms ), iAtomOut ) != end( group.atoms );
            if( !b_in_group ) return; // group cannot move if has rigid bond to something else
        }
        for( const size_t& iAtomOut : this->atoms[ iAtomIn ].bonded_atoms ) {
            bool b_in_group = find( begin( group.atoms ), end( group.atoms ), iAtomOut ) != end( group.atoms );
            if( b_in_group ) continue; 
            const Atom& atomIn  = this->atoms[ iAtomIn ];
            const Atom& atomOut = this->atoms[ iAtomOut ];
            if( !isWithinFlexibleBondNeighborhood( atomIn.x + dx, atomIn.y + dy, atomOut.x, atomOut.y ) ) {
                can_move = false;
                break;
            }
        }
    }
    if( !can_move ) return;
    // overlap test. 
    // simple implementation for now: remove from grid and try to place in the new position, else replace
    for( const auto& iAtom : group.atoms ) {
        const Atom &atom = this->atoms[ iAtom ];
        this->grid[ atom.x ][ atom.y ].has_atom = false;
    }
    bool all_ok = true;
    for( const auto& iAtom : group.atoms ) {
        const Atom &atom = this->atoms[ iAtom ];
        int tx = atom.x + dx;
        int ty = atom.y + dy;
        if( isOffGrid( tx, ty ) || this->grid[ tx ][ ty ].has_atom ) {
            all_ok = false;
            break;
        }
    }
    if( !all_ok ) {
        dx = dy = 0;
    }
    for( const auto& iAtom : group.atoms ) {
        Atom &atom = this->atoms[ iAtom ];
        atom.x += dx;
        atom.y += dy;
        this->grid[ atom.x ][ atom.y ].has_atom = true;
        this->grid[ atom.x ][ atom.y ].iAtom = iAtom;
    }
}

//----------------------------------------------------------------------------

