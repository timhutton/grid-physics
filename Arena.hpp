// STL:
#include <vector>

// Arena is a square grid world containing atoms
class Arena {

	public:
        
        Arena( int x, int y );

        // public typedefs
        struct Atom { int x, y; int type; std::vector<size_t> bonded_atoms; };
        enum BondType { Moore, vonNeumann };

		size_t addAtom( int x, int y, int type );
		void makeBond( size_t a, size_t b, BondType bond_type );
        void update();

        // accessors
        bool isOffGrid( int x, int y ) const;
        bool hasAtom( int x, int y ) const;
        int getArenaWidth() const { return this->X; }
        int getArenaHeight() const { return this->Y; }
        size_t getNumberOfAtoms() const { return this->atoms.size(); }
        Atom getAtom( size_t i ) const { return this->atoms[i]; }
        size_t getNumberOfGroups() const { return this->groups.size(); }
	
	private:

        // typedefs
        struct Group { std::vector<size_t> atoms; };
        struct Slot { bool has_atom; size_t iAtom; Slot() : has_atom( false ) {} };
        enum MovementType { JustAtoms      // atoms can move individually
                          , AllGroups      // all subgraphs of atoms can move individually
                          , MPEGSpace      // space itself moves around in large blocks
                          , MPEGMolecules  // molecules are divided spatially into movement blocks on the fly
                          };

        // private variables
        const int                         X;
        const int                         Y;
		std::vector<Group>                groups;
		std::vector<Atom>                 atoms;
        std::vector<std::vector<Slot>>    grid;
        const MovementType                movement_type;

        // private functions
        void addAllGroupsForNewBond( size_t a, size_t b );
        void removeGroupsWithOneButNotTheOther( size_t a, size_t b );
        bool moveGroupIfPossible( const Group& group, int dx, int dy );
        bool moveBlockIfPossible( int x, int y, int w, int h, int dx, int dy );
        void doChemistry();

        // useful constant values and functions
        static const int vNx[4];
        static const int vNy[4];
        static bool isWithinFlexibleBondNeighborhood( int x1, int y1, int x2, int y2 );
        static int getRandIntInclusive( int a, int b );
};
