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

size_t Arena::addAtom( int x, int y ) {

    if( isOffGrid(x,y ) )
		throw out_of_range("Atom not on grid");

	Slot& slot = this->grid[x][y];
	if( slot.has_atom )
		throw exception("Grid already contains an atom at that position");

	Atom a;
	a.x = x;
	a.y = y;
	this->atoms.push_back( a );
	size_t iAtom = this->atoms.size()-1;

    slot.has_atom = true;

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

	Bond bond;
	bond.a = a;
	bond.b = b;
	bond.type = type;
	this->bonds.push_back( bond );
	size_t iBond = this->bonds.size()-1;

	addFlexibleBond( a, b );
    if( type == BondType::vonNeumann ) {
		// remove any group that contains exactly one of a and b (can't have one without the other if a rigid bond)
		this->groups.erase( remove_if( begin( this->groups ), end( this->groups ), 
            GroupHasOneButNotTheOther( a, b ) ), end( this->groups ) );
    }
}

//----------------------------------------------------------------------------

void Arena::addFlexibleBond( size_t a, size_t b ) {
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
}

//----------------------------------------------------------------------------

bool Arena::isWithinFlexibleBondNeighborhood( int x1, int y1, int x2, int y2 ) {
    return abs( x1 - x2 ) <= 1 && abs( y1 - y2 ) <= 1;
}

//----------------------------------------------------------------------------

void Arena::update() {
    // attempt to move every group
    for( const auto& group : this->groups ) {
        int iMove = randInt( 4 );
        int dx = vNx[ iMove ];
        int dy = vNy[ iMove ];
        // first test: would this move stretch any bond too far?
        bool can_move = true;
        for( const auto& bond : this->bonds ) {
            bool a_in_group = find( begin( group.atoms ), end( group.atoms ), bond.a ) != end( group.atoms );
            bool b_in_group = find( begin( group.atoms ), end( group.atoms ), bond.b ) != end( group.atoms );
            if( ( a_in_group && b_in_group ) || (!a_in_group && !b_in_group ) )
                continue;
            size_t iAtomIn, iAtomOut;
            if( a_in_group ) { iAtomIn = bond.a; iAtomOut = bond.b; } else { iAtomIn = bond.b; iAtomOut = bond.a; }
            const Atom& atomIn  = this->atoms[ iAtomIn ];
            const Atom& atomOut = this->atoms[ iAtomOut ];
            if( !isWithinFlexibleBondNeighborhood( atomIn.x + dx, atomIn.y + dy, atomOut.x, atomOut.y ) ) {
                can_move = false;
                break;
            }
        }
        if( !can_move ) continue;
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
        }
    }
}

//----------------------------------------------------------------------------

