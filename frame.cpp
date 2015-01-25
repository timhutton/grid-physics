// local:
#include "frame.hpp"

// wxWidgets:
#include <wx/dcbuffer.h>

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
       , arena( 50, 50 )
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

    // an 8-cell loop with some rigid sections
    {
        size_t a = arena.addAtom( 1, 1 );
        size_t b = arena.addAtom( 2, 1 );
        size_t c = arena.addAtom( 2, 2 );
        size_t d = arena.addAtom( 1, 2 );
        size_t e = arena.addAtom( 1, 3 );
        size_t f = arena.addAtom( 0, 3 );
        size_t g = arena.addAtom( 0, 2 );
        size_t h = arena.addAtom( 0, 1 );
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
    {
        size_t a = arena.addAtom( 10, 10 );
        size_t b = arena.addAtom( 11, 10 );
        size_t c = arena.addAtom( 12, 10 );
        size_t d = arena.addAtom( 12, 11 );
        size_t e = arena.addAtom( 11, 11 );
        size_t f = arena.addAtom( 10, 11 );
        size_t g = arena.addAtom( 10, 12 );
        size_t h = arena.addAtom( 11, 12 );
        size_t i = arena.addAtom( 12, 12 );
        size_t j = arena.addAtom( 9, 9 );
        size_t k = arena.addAtom( 8, 8 );
        size_t l = arena.addAtom( 7, 7 );
        size_t m = arena.addAtom( 11, 9 );
        size_t n = arena.addAtom( 12, 8 );
        size_t o = arena.addAtom( 13, 7 );
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
    
    for( int i = 0; i < 100; ++i ) {
        int x = rand() % this->arena.getArenaWidth();
        int y = rand() % this->arena.getArenaHeight();
        if( !this->arena.hasAtom( x, y ) )
            this->arena.addAtom( x, y );
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

    pGC->SetPen(*wxMEDIUM_GREY_PEN);
    pGC->SetBrush(*wxLIGHT_GREY_BRUSH);
    for( size_t iAtom = 0; iAtom < this->arena.getNumberOfAtoms(); ++iAtom ) {
        Arena::Atom a = this->arena.getAtom( iAtom );
        pGC->DrawRectangle( a.x * scale, a.y * scale, scale, scale );
    }

    for( size_t iBond = 0; iBond < this->arena.getNumberOfBonds(); ++iBond) {
        Arena::Bond bond = this->arena.getBond( iBond );
        Arena::Atom a,b;
        a = this->arena.getAtom( bond.a );
        b = this->arena.getAtom( bond.b );
        switch( bond.type ) {
            case Arena::BondType::vonNeumann:  pGC->SetPen(*wxBLUE_PEN);       break;
            case Arena::BondType::Moore:       pGC->SetPen(*wxMEDIUM_GREY_PEN); break;
        }
        pGC->StrokeLine( ( a.x + 0.5 ) * scale, ( a.y + 0.5 ) * scale, ( b.x + 0.5 ) * scale, ( b.y + 0.5 ) * scale );
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

    this->arena.update();
    this->iterations++;
    if( iterations % this->render_every == 0 )
        this->Refresh( false );
    event.RequestMore(); // render continuously, not only once on idle
}

//-------------------------------------------------------------------------------------

void MyFrame::OnStep(wxCommandEvent& WXUNUSED(event)) {
    this->arena.update();
    this->iterations++;
    this->render_every = 0;
    this->Refresh(false);
}
