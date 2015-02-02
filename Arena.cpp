// local:
#include "Arena.hpp"

// stdlib
#include <limits.h>

// STL:
#include <stdexcept>
#include <algorithm>
using namespace std;

//----------------------------------------------------------------------------

Arena::Arena(int x, int y)
    : X( x )
	, Y( y )
    , movement_method( MovementMethod::MPEGMolecules )
    , movement_neighborhood( Neighborhood::vonNeumann ) // currently only vonNeumann supported
    , chemical_neighborhood( Neighborhood::vonNeumann )
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

size_t Arena::addAtom( int x, int y, int type ) {

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

void Arena::makeBond( size_t a, size_t b, Neighborhood range ) {
	if( a < 0 || a >= atoms.size() || b < 0 || b >= atoms.size() )
		throw out_of_range("Invalid atom index");
	if( a == b )
		throw invalid_argument("Cannot bond atom to itself");
    if( !isWithinNeighborhood( range, this->atoms[a].x, this->atoms[a].y, this->atoms[b].x, this->atoms[b].y ) )
        throw invalid_argument("Atoms are too far apart to be bonded");
    if( hasBond( a, b ) )
        throw invalid_argument("Atoms are already bonded");

    Bond ab = { b, range };
    this->atoms[ a ].bonds.push_back( ab );
    Bond ba = { a, range };
    this->atoms[ b ].bonds.push_back( ba );

    switch( this->movement_method ) {
        case JustAtoms:
            // here we can never move atoms with von Neumann bonds so
            // we can remove any groups with them in
            if( range == Neighborhood::vonNeumann )
                removeGroupsWithOneButNotTheOther( a, b );
            break;
        case AllGroups:
            //  we need to find all the subgraphs created by this new bond
	        addAllGroupsForNewBond( a, b );
            // if a and b are rigidly bonded then we don't need to check their groups separately
            if( range == Neighborhood::vonNeumann )
                removeGroupsWithOneButNotTheOther( a, b );
            break;
        case MPEGSpace: 
            // we don't use groups for this method
            break;
        case MPEGMolecules: 
            // here each molecule is a single group
            combineGroupsInvolvingTheseIntoOne( a, b );
            break;
    }
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

bool Arena::isWithinNeighborhood( Neighborhood type, int x1, int y1, int x2, int y2 ) {
    const int r2 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
    switch( type ) {
        case vonNeumann:  return r2 <= 1;
        case Moore:       return r2 <= 2;
        case vonNeumann2: return r2 <= 4;
        case knight:      return r2 <= 5;
        case Moore2:      return r2 <= 8;
    }
    throw out_of_range("unexpected enum");
}

//----------------------------------------------------------------------------

void Arena::getRandomMove( Neighborhood nhood, int &dx, int &dy ) {
    switch( nhood ) {
        case Neighborhood::vonNeumann: {
            const int vNx[4] = {  0,  1,  0, -1 }; // clockwise from North
            const int vNy[4] = { -1,  0,  1,  0 };
            const int iMove = getRandIntInclusive(0,3);
            dx = vNx[ iMove ];
            dy = vNy[ iMove ];
            break;
        }
        case Neighborhood::Moore: {
            const int Mx[8] = {  0,  1,  1,  1,  0, -1, -1, -1 }; // clockwise from North
            const int My[8] = { -1, -1,  0,  1,  1,  1,  0, -1 };
            const int iMove = getRandIntInclusive(0,7);
            dx = Mx[ iMove ];
            dy = My[ iMove ];
            break;
        }
        default: throw out_of_range("Unsupported neighborhood");
    }
}

//----------------------------------------------------------------------------

void Arena::update() {
    switch( this->movement_method ) {
        case JustAtoms:
        case AllGroups:
            // attempt to move every group
            for( const auto& group : this->groups ) {
                int dx, dy;
                getRandomMove( this->movement_neighborhood, dx, dy );
                moveGroupIfPossible( group, dx, dy );
            }
            break;
        case MPEGSpace: {
            for( int i = 0; i < 10; ++i ) {
                // attempt to move a block
                int x = getRandIntInclusive( 0, this->X-1 );
                int y = getRandIntInclusive( 0, this->Y-1 );
                int w = getRandIntInclusive( 1, this->X-x );
                int h = getRandIntInclusive( 1, this->Y-y );
                int dx, dy;
                getRandomMove( this->movement_neighborhood, dx, dy );
                moveBlockIfPossible( x, y, w, h, dx, dy );
            }
            break;
        }
        case MPEGMolecules:
            // attempt to move every group
            for( const auto& group : this->groups ) {
                moveBlocksInGroup( group );
            }
            break;
    }

    // find chemical reactions
    doChemistry();
}

//----------------------------------------------------------------------------

void Arena::doChemistry() {
    for( int x = 0; x < this->X; ++x ) {
        for( int y = 0; y < this->Y; ++y ) {
            if( !this->grid[x][y].has_atom ) continue;
            int dx, dy;
            getRandomMove( this->chemical_neighborhood, dx, dy );
            int tx = x + dx;
            int ty = y + dy;
            if( isOffGrid( tx, ty ) || !this->grid[tx][ty].has_atom ) continue;
            size_t iAtomA = this->grid[x][y].iAtom;
            size_t iAtomB = this->grid[tx][ty].iAtom;
            Atom& a = this->atoms[ iAtomA ];
            Atom& b = this->atoms[ iAtomB ];
            if( !hasBond( iAtomA, iAtomB ) && a.type == b.type && a.bonds.size() + b.bonds.size() < 2 ) {
                Neighborhood bond_range = Neighborhood::Moore;
                makeBond( iAtomA, iAtomB, bond_range );
            }
        }
    }
}

//----------------------------------------------------------------------------

bool Arena::moveGroupIfPossible( const Group& group, int dx, int dy ) {
    // first test: would this move stretch any bond too far?
    bool can_move = true;
    for( const size_t& iAtomIn : group.atoms ) {
        for( const Bond& bond : this->atoms[ iAtomIn ].bonds ) {
            const size_t& iAtomOut = bond.iAtom;
            bool b_in_group = find( begin( group.atoms ), end( group.atoms ), iAtomOut ) != end( group.atoms );
            if( b_in_group ) continue; 
            const Atom& atomIn  = this->atoms[ iAtomIn ];
            const Atom& atomOut = this->atoms[ iAtomOut ];
            if( !isWithinNeighborhood( bond.range, atomIn.x + dx, atomIn.y + dy, atomOut.x, atomOut.y ) ) {
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
            for( const Bond& bond : a.bonds ) {
                const size_t iAtomB = bond.iAtom;
                const Atom& b = this->atoms[ iAtomB ];
                if( b.x >= left && b.x <= right && b.y >= top && b.y <= bottom )
                    continue; // atom B is also within the block
                if( !isWithinNeighborhood( bond.range, sx + dx, sy + dy, b.x, b.y ) )
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

void Arena::moveBlocksInGroup( const Group& group ) {
    // get the bounding box
    int bb[4] = { INT_MAX, -INT_MAX, INT_MAX, -INT_MAX };
    for( const size_t& iAtom : group.atoms ) {
        const Atom& a = this->atoms[ iAtom ];
        bb[0] = min( bb[0], a.x );
        bb[1] = max( bb[1], a.x );
        bb[2] = min( bb[2], a.y );
        bb[3] = max( bb[3], a.y );
    }
    // let the whole block have a go at moving
    int dx, dy;
    getRandomMove( this->movement_neighborhood, dx, dy );
    bool moved = moveMembersOfGroupInBlockIfPossible( group, bb[0], bb[2], bb[1]-bb[0]+1, bb[3]-bb[2]+1, dx, dy );
    if( moved ) { bb[0]+=dx; bb[1]+=dx; bb[2]+=dy; bb[3]+=dy; }
    // also try moving some rectangles within it
    // (some bits might move outside the bounding box but that's OK)
    const int num_tries = ( bb[1] - bb[0] + 1) * ( bb[3] - bb[2] + 1 );
    for( int iTry = 0; iTry < num_tries; ++iTry ) {
        int w = getRandIntInclusive( 1, bb[1] - bb[0] + 1 );
        int h = getRandIntInclusive( 1, bb[3] - bb[2] + 1 );
        int x = getRandIntInclusive( bb[0], bb[1] - w + 1 );
        int y = getRandIntInclusive( bb[2], bb[3] - h + 1 );
        getRandomMove( this->movement_neighborhood, dx, dy );
        moveMembersOfGroupInBlockIfPossible( group, x, y, w, h, dx, dy );
    }
}

//----------------------------------------------------------------------------

bool Arena::moveMembersOfGroupInBlockIfPossible( const Group& group, int x, int y, int w, int h, int dx, int dy ) {
    // collect the atoms in this block that we want to move
    vector<size_t> movers;
    for( int sx = x; sx < x+w; ++sx ) {
        for( int sy = y; sy < y+h; ++sy ) {
            if( isOffGrid( sx, sy ) || !this->grid[sx][sy].has_atom )
                continue; // not an atom here
            const size_t iAtom = this->grid[sx][sy].iAtom;
            if( find( group.atoms.begin(), group.atoms.end(), iAtom ) == group.atoms.end() )
                continue; // not one of our group's atoms
            int tx = sx + dx;
            int ty = sy + dy;
            if( isOffGrid( tx, ty ) )
                return false; // can't move off-grid
            movers.push_back( iAtom );
        }
    }
    // bond check
    for( const size_t& iAtom : movers ) {
        const Atom& a = this->atoms[ iAtom ];
        for( const Bond& bond : a.bonds ) {
            const size_t iAtomB = bond.iAtom;
            if( find( movers.begin(), movers.end(), iAtomB ) != movers.end() )
                continue; // no problem, since B is also part of the moving set
            const Atom& b = this->atoms[ iAtomB ];
            if( !isWithinNeighborhood( bond.range, a.x + dx, a.y + dy, b.x, b.y ) )
                return false; // would over-stretch this bond
        }
    }
    // overlap check: 
    // simple implementation for now: remove from grid and try to place in the new position, else replace
    for( const size_t& iAtom : movers ) {
        const Atom &atom = this->atoms[ iAtom ];
        this->grid[ atom.x ][ atom.y ].has_atom = false;
    }
    bool all_ok = true;
    for( const size_t& iAtom : movers ) {
        const Atom &a = this->atoms[ iAtom ];
        int tx = a.x + dx;
        int ty = a.y + dy;
        if( this->grid[ tx ][ ty ].has_atom ) {
            all_ok = false;
            break;
        }
    }
    if( !all_ok ) {
        dx = dy = 0;
    }
    for( const size_t& iAtom : movers ) {
        Atom &a = this->atoms[ iAtom ];
        a.x += dx;
        a.y += dy;
        this->grid[ a.x ][ a.y ].has_atom = true;
        this->grid[ a.x ][ a.y ].iAtom = iAtom;
    }
    return all_ok;
}
                                
//----------------------------------------------------------------------------

int Arena::getRandIntInclusive( int a, int b )
{
    return a + rand() % ( b - a + 1 );
}

//----------------------------------------------------------------------------

void Arena::combineGroupsInvolvingTheseIntoOne( size_t a, size_t b ) {
    // find every group involving a or b
    vector<size_t> groups_to_be_merged;
	for( size_t iGroup = 0; iGroup < this->groups.size(); ++iGroup ) {
        const Group& g = this->groups[ iGroup ];
		if( find( begin( g.atoms ), end( g.atoms ), a ) != end( g.atoms ) ||
                find( begin( g.atoms ), end( g.atoms ), b ) != end( g.atoms ) )
            groups_to_be_merged.push_back( iGroup );
    }
    if( groups_to_be_merged.size() < 2 )
        return; // nothing to do

	Group& g = this->groups[ groups_to_be_merged.front() ];
    for( size_t iiGroup = 1; iiGroup < groups_to_be_merged.size(); ++iiGroup ) {
        const Group& gb = this->groups[ groups_to_be_merged[ iiGroup ] ];
        vector<size_t> merged_atoms( g.atoms.size() + gb.atoms.size() );
        const auto& end = set_union( g.atoms.begin(), g.atoms.end(), gb.atoms.begin(), gb.atoms.end(), merged_atoms.begin() );
        merged_atoms.resize( end - merged_atoms.begin() );
        g.atoms.assign( merged_atoms.begin(), merged_atoms.end() );
    }
    for( size_t iiGroup = groups_to_be_merged.size()-1; iiGroup >= 1; --iiGroup ) {
        this->groups.erase( this->groups.begin() + groups_to_be_merged[ iiGroup  ] );
    }
}

//----------------------------------------------------------------------------

bool Arena::hasBond( size_t a, size_t b ) const {
    for( const Bond& bond : this->atoms[ a ].bonds ) {
        if( bond.iAtom == b )
            return true;
    }
    return false;
}

//----------------------------------------------------------------------------
