// local:
#include "Arena.hpp"

// STL:
#include <stdexcept>
#include <algorithm>
using namespace std;

//----------------------------------------------------------------------------

const int Arena::vNx[4] = { 0, 1, 0, -1 }; // NESW
const int Arena::vNy[4] = { -1, 0, 1, 0 };

//----------------------------------------------------------------------------

Arena::Arena(int x, int y)
    : X( x )
	, Y( y )
    , movement_type( MovementType::JustAtoms )
{
	this->grid = vector<vector<Slot>>( X, vector<Slot>( Y ) );
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

    this->atoms[ a ].bonded_atoms.push_back( b );
    this->atoms[ b ].bonded_atoms.push_back( a );

    if( this->movement_type == MovementType::AllGroups )
	    addAllGroupsForFlexibleBond( a, b );

    if( type == BondType::vonNeumann )
        removeGroupsWithOneButNotTheOther( a, b );
}

//----------------------------------------------------------------------------

void Arena::addAllGroupsForFlexibleBond( size_t a, size_t b ) {
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

void Arena::removeGroupsWithOneButNotTheOther( size_t a, size_t b ) {

    class GroupHasOneButNotTheOther {
        public:
            GroupHasOneButNotTheOther( size_t a, size_t b ) : a(a), b(b) {}
            bool operator() (const Group& g) const
            { 
                int n = 0;
                if( find( begin(g.atoms), end(g.atoms), a) != end(g.atoms) ) n++;
                if( find( begin(g.atoms), end(g.atoms), b) != end(g.atoms) ) n++;
                return n == 1;
            }
        private:
            size_t a,b;
    };

    this->groups.erase( remove_if( begin( this->groups ), end( this->groups ), 
        GroupHasOneButNotTheOther( a, b ) ), end( this->groups ) );
}

//----------------------------------------------------------------------------

bool Arena::isWithinFlexibleBondNeighborhood( int x1, int y1, int x2, int y2 ) {
    return abs( x1 - x2 ) <= 1 && abs( y1 - y2 ) <= 1;
}

//----------------------------------------------------------------------------

void Arena::update() {
    // attempt to move every group
    for( const auto& group : this->groups ) {
        moveGroupRandomly( group );
    }
    // find chemical reactions
    doChemistry();
}

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
    int iMove = getRandInt( 4 );
    int dx = vNx[ iMove ];
    int dy = vNy[ iMove ];
    // first test: would this move stretch any bond too far?
    bool can_move = true;
    for( const size_t& iAtomIn : group.atoms ) {
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

int Arena::getRandInt( int n )  // [0,n)
{
    return rand() % n; // use something more uniform if needed
}

//----------------------------------------------------------------------------
