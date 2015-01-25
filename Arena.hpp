// STL:
#include <vector>
using namespace std;

// Arena is a square grid world containing atoms
class Arena {

	public:
        
        struct Atom { int x, y; vector<size_t> bonded_atoms; };
        enum BondType { Moore, vonNeumann };

        Arena( int x, int y );

        bool isOffGrid( int x, int y ) const;
        bool hasAtom( int x, int y ) const;
		size_t addAtom( int x, int y );
		void makeBond( size_t a, size_t b, BondType bond_type );
        void update();

        int getArenaWidth() const { return this->X; }
        int getArenaHeight() const { return this->Y; }

        size_t getNumberOfAtoms() const { return this->atoms.size(); }
        Atom getAtom( size_t i ) const { return this->atoms[i]; }
	
        size_t getNumberOfGroups() const { return this->groups.size(); }
	
	private:

        struct Group { vector<size_t> atoms; };
        struct Slot { bool has_atom; size_t iAtom; Slot() : has_atom( false ) {} };

        int X;
        int Y;
		vector<Group>        groups;
		vector<Atom>         atoms;
        vector<vector<Slot>> grid;

        void addFlexibleBond( size_t a, size_t b );
        static bool isWithinFlexibleBondNeighborhood( int x1, int y1, int x2, int y2 );
        void moveGroupRandomly( const Group& group );
        void doChemistry();

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

        static const int vNx[4];
        static const int vNy[4];
};
