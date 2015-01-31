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
    , movement_type( MovementType::AllGroups )
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
	    addAllGroupsForNewBond( a, b );

    if( type == BondType::vonNeumann )
        removeGroupsWithOneButNotTheOther( a, b );
}

//----------------------------------------------------------------------------

void Arena::addAllGroupsForNewBond( size_t a, size_t b ) {
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
    switch( this->movement_type ) {
        case JustAtoms:
        case AllGroups:
            // attempt to move every group
            for( const auto& group : this->groups ) {
                int iMove = getRandIntInclusive( 0, 3 );
                moveGroupIfPossible( group, vNx[ iMove ], vNy[ iMove ] );
            }
            break;
        case MPEGSpace: {
            // attempt to move a block
            int x = getRandIntInclusive( 0, this->X-1 );
            int y = getRandIntInclusive( 0, this->Y-1 );
            int w = getRandIntInclusive( 1, this->X-x );
            int h = getRandIntInclusive( 1, this->Y-y );
            int iMove = getRandIntInclusive( 0, 3 );
            moveBlockIfPossible( x, y, w, h, vNx[ iMove ], vNy[ iMove ] );
            break;
        }
        case MPEGMolecules:
            throw exception("Not yet implemented");
            break;
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

bool Arena::moveGroupIfPossible( const Group& group, int dx, int dy ) {
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
    if( !can_move ) return false;
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
    return all_ok;
}

//----------------------------------------------------------------------------

bool Arena::moveBlockIfPossible( int x, int y, int w, int h, int dx, int dy ) {
    const int left = x;
    const int right = x + w - 1;
    const int top = y;
    const int bottom = y + h - 1;
    if( isOffGrid( left, top ) || isOffGrid( right, bottom ) )
        throw out_of_range("Attempt to move block that is not wholy on the grid");
    // overlap test along front edge:
    int x1, y1, x2, y2;
    if( dx == 1 )       { x1 = x2 = right;  y1 = top;  y2 = bottom; }
    else if( dx == -1 ) { x1 = x2 = left;   y1 = top;  y2 = bottom; }
    else if( dy == 1 )  { y1 = y2 = bottom; x1 = left; x2 = right;  }
    else if( dy == -1 ) { y1 = y2 = top;    x1 = left; x2 = right;  }
    for( int sx = x1; sx <= x2; ++sx ) {
        for( int sy = y1; sy <= y2; ++sy ) {
            if( !this->grid[sx][sy].has_atom )
                continue;
            int tx = sx + dx;
            int ty = sy + dy;
            if( isOffGrid(tx,ty) || this->grid[tx][ty].has_atom )
                return false;
        }
    }
    // bond test along back and side edges:
    // TODO: improve test to reject leading edge
    // TODO: find better way of traversing non-leading edges
    for( int sx = left; sx <= right; ++sx ) {
        for( int sy = top; sy <= bottom; ++sy ) {
            if( sx > left && sx < right && sy > top && sy < bottom )
                continue; // not on the edge of the block
            if( !this->grid[sx][sy].has_atom )
                continue;
            const Atom& a = this->atoms[ this->grid[sx][sy].iAtom ];
            for( const size_t iAtomB : a.bonded_atoms ) {
                const Atom& b = this->atoms[ iAtomB ];
                if( b.x >= left && b.x <= right && b.y >= top && b.y <= bottom )
                    continue; // atom B is also within the block
                if( !isWithinFlexibleBondNeighborhood( sx + dx, sy + dy, b.x, b.y ) )
                    return false; // would over-stretch this bond
            }
        }
    }
    // move the block
    vector<size_t> movers;
    for( int sx = left; sx <= right; ++sx ) {
        for( int sy = top; sy <= bottom; ++sy ) {
            if( !this->grid[sx][sy].has_atom )
                continue;
            movers.push_back( this->grid[sx][sy].iAtom );
            this->grid[sx][sy].has_atom = false;
        }
    }
    for( const size_t& iAtom : movers ) {
        Atom& a = this->atoms[ iAtom ];
        a.x += dx;
        a.y += dy;
        if( isOffGrid( a.x, a.y ) )
            throw exception("internal error");
        Slot& slot = this->grid[ a.x ][ a.y ];
        slot.has_atom = true;
        slot.iAtom = iAtom;
    }
    return true;
}

//----------------------------------------------------------------------------

int Arena::getRandIntInclusive( int a, int b )
{
    return a + rand() % (b-a+1); // use something more uniform if needed
}

//----------------------------------------------------------------------------
