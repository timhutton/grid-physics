// STL:
#include <vector>

// Arena is a square grid world containing atoms
class Arena {

	public:
        
        Arena( int x, int y );

        // public typedefs                  // as squared Euclidean distance r2:
        enum Neighborhood { vonNeumann      // r2 <= 1
                          , Moore           // r2 <= 2
                          , vonNeumann2     // r2 <= 4
                          , knight          // r2 <= 5
                          , Moore2          // r2 <= 8
                          };
        struct Bond { size_t iAtom; Neighborhood range; };
        struct Atom { int x, y; int type; std::vector<Bond> bonds; };

		size_t addAtom( int x, int y, int type );
		void makeBond( size_t a, size_t b, Neighborhood range );
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
        enum MovementMethod { JustAtoms      // atoms can move individually
                            , AllGroups      // all subgraphs of atoms can move individually
                            , MPEGSpace      // space itself moves around in large blocks
                            , MPEGMolecules  // molecules are divided spatially into movement blocks on the fly
                            };

        // private variables
        const int                         X;
        const int                         Y;
		std::vector<Atom>                 atoms;
        std::vector<std::vector<Slot>>    grid;
		std::vector<Group>                groups;
        const MovementMethod              movement_method;
        const Neighborhood                movement_neighborhood;
        const Neighborhood                chemical_neighborhood;

        // private functions
        void addAllGroupsForNewBond( size_t a, size_t b );
        void removeGroupsWithOneButNotTheOther( size_t a, size_t b );
        void combineGroupsInvolvingTheseIntoOne( size_t a, size_t b );
        bool moveGroupIfPossible( const Group& group, int dx, int dy );
        bool moveBlockIfPossible( int x, int y, int w, int h, int dx, int dy );
        void moveBlocksInGroup( const Group& group );
        void moveBlocksInGroup( const Group& group, int x, int y, int w, int h );
        bool moveMembersOfGroupInBlockIfPossible( const Group& group, int x, int y, int w, int h, int dx, int dy  );
        void doChemistry();
        bool hasBond( size_t a, size_t b ) const;
        int getRandomMove() const;

        // useful functions
        static bool isWithinNeighborhood( Neighborhood type, int x1, int y1, int x2, int y2 );
        static int getRandIntInclusive( int a, int b );
        static void getRandomMove( Neighborhood nhood, int& dx, int& dy );
};
