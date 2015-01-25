// local:
#include "frame.hpp"

// wxWidgets:
#include <wx/dcbuffer.h>

wxBEGIN_EVENT_TABLE(MyFrame, wxFrame)
    EVT_MENU(wxID_EXIT,  MyFrame::OnQuit)
    EVT_MENU(wxID_ABOUT, MyFrame::OnAbout)
    EVT_PAINT(MyFrame::OnPaint)
    EVT_SIZE(MyFrame::OnSize)
    EVT_IDLE(MyFrame::OnIdle)
wxEND_EVENT_TABLE()

//-------------------------------------------------------------------------------------

MyFrame::MyFrame(const wxString& title)
       : wxFrame(NULL, wxID_ANY, title)
       , arena( 50, 50 )
       , iterations( 0 )
       , render_every( 100 )
{
    SetIcon(wxICON(sample));

#if wxUSE_MENUS
    wxMenu *fileMenu = new wxMenu;
    fileMenu->Append(wxID_EXIT, "E&xit\tAlt-X", "Quit this program");

    wxMenu *helpMenu = new wxMenu;
    helpMenu->Append(wxID_ABOUT, "&About\tF1", "Show about dialog");

    wxMenuBar *menuBar = new wxMenuBar();
    menuBar->Append(fileMenu, "&File");
    menuBar->Append(helpMenu, "&Help");

    SetMenuBar(menuBar);
#endif // wxUSE_MENUS

    SetBackgroundStyle( wxBackgroundStyle::wxBG_STYLE_PAINT );

    seed();
}

//-------------------------------------------------------------------------------------

void MyFrame::seed() {
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
    
    for( int i = 0; i < 20; ++i ) {
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
    pGC->DrawText( wxString::Format("%d",this->iterations), 10, 10 );
}

//-------------------------------------------------------------------------------------

void MyFrame::OnSize(wxSizeEvent& event) {
    this->Refresh( false );
}

//-------------------------------------------------------------------------------------

void MyFrame::OnIdle(wxIdleEvent& event) {
    this->arena.update();
    this->iterations++;
    if( iterations % this->render_every == 0 )
        this->Refresh( false );
    event.RequestMore(); // render continuously, not only once on idle
}

//-------------------------------------------------------------------------------------

