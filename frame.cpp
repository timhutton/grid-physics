// local:
#include "frame.hpp"

// wxWidgets:
#include <wx/dcbuffer.h>

// STL:
using namespace std;

namespace ID
{
    enum
    {
        // we can use IDs higher than this for our own purposes
        Dummy = wxID_HIGHEST+1,

        SpeedStop,
        SpeedSlowest,
        SpeedMedium,
        SpeedFast,
        Step,
    };
};

wxBEGIN_EVENT_TABLE(MyFrame, wxFrame)
    EVT_MENU(wxID_EXIT,  MyFrame::OnQuit)
    EVT_MENU(wxID_ABOUT, MyFrame::OnAbout)
    EVT_MENU(ID::SpeedStop, MyFrame::OnSpeedStop)
    EVT_MENU(ID::SpeedSlowest, MyFrame::OnSpeedSlowest)
    EVT_MENU(ID::SpeedMedium, MyFrame::OnSpeedMedium)
    EVT_MENU(ID::SpeedFast, MyFrame::OnSpeedFast)
    EVT_MENU(ID::Step, MyFrame::OnStep)
    EVT_PAINT(MyFrame::OnPaint)
    EVT_SIZE(MyFrame::OnSize)
    EVT_IDLE(MyFrame::OnIdle)
wxEND_EVENT_TABLE()

//-------------------------------------------------------------------------------------

MyFrame::MyFrame(const wxString& title)
       : wxFrame(NULL, wxID_ANY, title, wxDefaultPosition, wxSize(900,700) )
       , arena( 80, 60 )
       , iterations( 0 )
       , render_every( 1 )
{
    SetIcon(wxICON(sample));

#if wxUSE_MENUS
    wxMenu *fileMenu = new wxMenu;
    fileMenu->Append(wxID_EXIT, "E&xit\tAlt-X", "Quit this program");

    wxMenu *actionMenu = new wxMenu;
    actionMenu->Append(ID::SpeedStop, "Stop\t0", "Stop running");
    actionMenu->Append(ID::SpeedSlowest, "Run at slow speed\t1", "Run at the slowest speed");
    actionMenu->Append(ID::SpeedMedium, "Run at medium speed\t2", "Run at a medium speed");
    actionMenu->Append(ID::SpeedFast, "Run at fast speed\t3", "Run at a fast speed");
    actionMenu->Append(ID::Step, "Step\tSPACE", "Advance forwards by a single step");

    wxMenu *helpMenu = new wxMenu;
    helpMenu->Append(wxID_ABOUT, "&About\tF1", "Show about dialog");

    wxMenuBar *menuBar = new wxMenuBar();
    menuBar->Append(fileMenu, "&File");
    menuBar->Append(actionMenu, "&Action");
    menuBar->Append(helpMenu, "&Help");

    SetMenuBar(menuBar);
#endif // wxUSE_MENUS

    SetBackgroundStyle( wxBackgroundStyle::wxBG_STYLE_PAINT );

    seed();
}

//-------------------------------------------------------------------------------------

void MyFrame::seed() {

    try {
        // an 8-cell loop with some rigid sections
        if( 1 ) {
            size_t a = arena.addAtom( 1, 1, 0 );
            size_t b = arena.addAtom( 2, 1, 0 );
            size_t c = arena.addAtom( 2, 2, 0 );
            size_t d = arena.addAtom( 1, 2, 0 );
            size_t e = arena.addAtom( 1, 3, 0 );
            size_t f = arena.addAtom( 0, 3, 0 );
            size_t g = arena.addAtom( 0, 2, 0 );
            size_t h = arena.addAtom( 0, 1, 0 );
            arena.makeBond( a, b, Arena::BondType::vonNeumann );
            arena.makeBond( b, c, Arena::BondType::vonNeumann );
            arena.makeBond( c, d, Arena::BondType::Moore );
            arena.makeBond( d, e, Arena::BondType::Moore );
            arena.makeBond( e, f, Arena::BondType::vonNeumann );
            arena.makeBond( f, g, Arena::BondType::Moore );
            arena.makeBond( g, h, Arena::BondType::vonNeumann );
            arena.makeBond( h, a, Arena::BondType::Moore );
        }
    
        // a box with flailing arms
        if( 1 ) {
            size_t a = arena.addAtom( 10, 10, 1 );
            size_t b = arena.addAtom( 11, 10, 1 );
            size_t c = arena.addAtom( 12, 10, 1 );
            size_t d = arena.addAtom( 12, 11, 1 );
            size_t e = arena.addAtom( 11, 11, 1 );
            size_t f = arena.addAtom( 10, 11, 1 );
            size_t g = arena.addAtom( 10, 12, 1 );
            size_t h = arena.addAtom( 11, 12, 1 );
            size_t i = arena.addAtom( 12, 12, 1 );
            size_t j = arena.addAtom( 9, 9, 1 );
            size_t k = arena.addAtom( 8, 8, 1 );
            size_t l = arena.addAtom( 7, 7, 1 );
            size_t m = arena.addAtom( 11, 9, 1 );
            size_t n = arena.addAtom( 12, 8, 1 );
            size_t o = arena.addAtom( 13, 7, 1 );
            arena.makeBond( a, b, Arena::BondType::vonNeumann );
            arena.makeBond( b, c, Arena::BondType::vonNeumann );
            arena.makeBond( c, d, Arena::BondType::vonNeumann );
            arena.makeBond( d, e, Arena::BondType::vonNeumann );
            arena.makeBond( e, f, Arena::BondType::vonNeumann );
            arena.makeBond( f, g, Arena::BondType::vonNeumann );
            arena.makeBond( g, h, Arena::BondType::vonNeumann );
            arena.makeBond( h, i, Arena::BondType::vonNeumann );
            arena.makeBond( a, j, Arena::BondType::Moore );
            arena.makeBond( j, k, Arena::BondType::Moore );
            arena.makeBond( k, l, Arena::BondType::Moore );
            arena.makeBond( c, m, Arena::BondType::Moore );
            arena.makeBond( m, n, Arena::BondType::Moore );
            arena.makeBond( n, o, Arena::BondType::Moore );
        }

        // a double-stranded molecule
        if( 1 ) {
            size_t a = arena.addAtom( 21, 21, 5 );
            size_t b = arena.addAtom( 22, 21, 3 );
            size_t c = arena.addAtom( 21, 22, 5 );
            size_t d = arena.addAtom( 22, 22, 3 );
            size_t e = arena.addAtom( 21, 23, 5 );
            size_t f = arena.addAtom( 22, 23, 3 );
            size_t g = arena.addAtom( 21, 24, 5 );
            size_t h = arena.addAtom( 22, 24, 3 );
            size_t i = arena.addAtom( 21, 25, 5 );
            size_t j = arena.addAtom( 22, 25, 3 );
            size_t k = arena.addAtom( 21, 26, 5 );
            size_t l = arena.addAtom( 22, 26, 3 );
            arena.makeBond( a, b, Arena::BondType::Moore );
            arena.makeBond( c, d, Arena::BondType::Moore );
            arena.makeBond( e, f, Arena::BondType::Moore );
            arena.makeBond( g, h, Arena::BondType::Moore );
            arena.makeBond( i, j, Arena::BondType::Moore );
            arena.makeBond( k, l, Arena::BondType::Moore );
            arena.makeBond( a, c, Arena::BondType::Moore );
            arena.makeBond( b, d, Arena::BondType::Moore );
            arena.makeBond( c, e, Arena::BondType::Moore );
            arena.makeBond( d, f, Arena::BondType::Moore );
            arena.makeBond( e, g, Arena::BondType::Moore );
            arena.makeBond( f, h, Arena::BondType::Moore );
            arena.makeBond( g, i, Arena::BondType::Moore );
            arena.makeBond( h, j, Arena::BondType::Moore );
            arena.makeBond( i, k, Arena::BondType::Moore );
            arena.makeBond( j, l, Arena::BondType::Moore );
        }

        if( 1 ) { 
            // a longer chain
            const int N = 10;
            size_t a = arena.addAtom( 31, 0, 2 );
            size_t b = arena.addAtom( 32, 0, 2 );
            arena.makeBond( a, b, Arena::BondType::Moore );
            for( int i = 1; i < N; ++i ) {
                size_t a2 = arena.addAtom( 31, i, 2 );
                size_t b2 = arena.addAtom( 32, i, 2 );
                arena.makeBond( a, a2, Arena::BondType::Moore );
                arena.makeBond( b, b2, Arena::BondType::Moore );
                a = a2;
                b = b2;
            }
        }
    
        if( 1 ) {
            // add some surrounding atoms
            for( int i = 0; i < 500; ++i ) {
                int x = rand() % this->arena.getArenaWidth();
                int y = rand() % this->arena.getArenaHeight();
                if( !this->arena.hasAtom( x, y ) )
                    this->arena.addAtom( x, y, rand() % 6 );
            }
        }
    }
    catch( exception& e ) {
        wxMessageBox( e.what() );
    }
}

//-------------------------------------------------------------------------------------

void MyFrame::OnQuit(wxCommandEvent& WXUNUSED(event))
{
    Close(true);
}

//-------------------------------------------------------------------------------------

void MyFrame::OnAbout(wxCommandEvent& WXUNUSED(event))
{
    wxMessageBox("Grid Physics",
                 "About Grid Physics",
                 wxOK | wxICON_INFORMATION,
                 this);
}

//-------------------------------------------------------------------------------------

void MyFrame::OnPaint(wxPaintEvent& WXUNUSED(event))
{
    int X = this->GetClientSize().x; 
    int Y = this->GetClientSize().y; 
    if( X < 0 || Y < 0 )
        return;

    wxAutoBufferedPaintDC dc( this );        
    wxGraphicsContext* pGC = wxGraphicsContext::Create( dc );
    if( pGC ) {
        draw( pGC, X, Y );
    }
}

//-------------------------------------------------------------------------------------

void MyFrame::draw( wxGraphicsContext *pGC, int X, int Y ) {
    int scale = min( (X-1) / this->arena.getArenaWidth(), (Y-1) / this->arena.getArenaHeight() );

    // blank the area
    pGC->SetPen(*wxLIGHT_GREY_PEN);
    pGC->SetBrush(*wxLIGHT_GREY_BRUSH);
    pGC->DrawRectangle( 0, 0, X, Y );

    // draw the arena
    pGC->SetPen(*wxBLACK_PEN);
    pGC->SetBrush(*wxWHITE_BRUSH);
    pGC->DrawRectangle( 0, 0, this->arena.getArenaWidth()*scale, this->arena.getArenaHeight()*scale );
    drawArena( pGC, scale );

    delete pGC;
}

//-------------------------------------------------------------------------------------

void MyFrame::drawArena( wxGraphicsContext* pGC, int scale ) {

    // draw atoms
    pGC->SetPen(*wxMEDIUM_GREY_PEN);
    pGC->SetBrush(*wxLIGHT_GREY_BRUSH);
    for( size_t iAtom = 0; iAtom < this->arena.getNumberOfAtoms(); ++iAtom ) {
        Arena::Atom a = this->arena.getAtom( iAtom );
        switch( a.type ) {
            default:
            case 0: pGC->SetBrush(*wxRED_BRUSH); break;
            case 1: pGC->SetBrush(*wxYELLOW_BRUSH); break;
            case 2: pGC->SetBrush(*wxCYAN_BRUSH); break;
            case 3: pGC->SetBrush(*wxLIGHT_GREY_BRUSH); break;
            case 4: pGC->SetBrush(*wxBLUE_BRUSH); break;
            case 5: pGC->SetBrush(*wxGREEN_BRUSH); break;
        }
        pGC->DrawRectangle( a.x * scale, a.y * scale, scale, scale );
    }

    // draw bonds
    wxPen thinBondPen(*wxBLACK,1);
    wxPen thickBondPen(*wxBLACK,2);
    for( size_t iAtom = 0; iAtom < this->arena.getNumberOfAtoms(); ++iAtom ) {
        Arena::Atom a = this->arena.getAtom( iAtom );
        for( const Arena::Bond& bond : a.bonds ) {
            const size_t iAtom2 = bond.iAtom;
            if( iAtom2 < iAtom ) continue; 
            Arena::Atom b = this->arena.getAtom( iAtom2 );
            switch( bond.type ) {
                case Arena::BondType::Moore:      pGC->SetPen(thinBondPen);  break;
                case Arena::BondType::vonNeumann: pGC->SetPen(thickBondPen); break;
            }
            pGC->StrokeLine( ( a.x + 0.5 ) * scale, ( a.y + 0.5 ) * scale, ( b.x + 0.5 ) * scale, ( b.y + 0.5 ) * scale );
        }
    }

    pGC->SetFont(*wxNORMAL_FONT,*wxBLACK);
    pGC->DrawText( wxString::Format("Its: %d",this->iterations), 10, 10 );
    pGC->DrawText( wxString::Format("Groups: %d",this->arena.getNumberOfGroups()), 10, 30 );
}

//-------------------------------------------------------------------------------------

void MyFrame::OnSize(wxSizeEvent& event) {
    this->Refresh( false );
}

//-------------------------------------------------------------------------------------

void MyFrame::OnIdle(wxIdleEvent& event) {
    if( render_every == 0 ) return;

    try {
        this->arena.update();
    }
    catch( exception& e ) {
        wxMessageBox( e.what() );
    }
    this->iterations++;
    if( iterations % this->render_every == 0 )
        this->Refresh( false );
    event.RequestMore(); // render continuously, not only once on idle
}

//-------------------------------------------------------------------------------------

void MyFrame::OnStep(wxCommandEvent& WXUNUSED(event)) {
    try {
        this->arena.update();
    }
    catch( exception& e ) {
        wxMessageBox( e.what() );
    }
    this->iterations++;
    this->render_every = 0;
    this->Refresh(false);
}
