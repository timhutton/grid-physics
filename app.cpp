// local:
#include "app.hpp"
#include "frame.hpp"

IMPLEMENT_APP(MyApp)

bool MyApp::OnInit()
{
    if ( !wxApp::OnInit() )
        return false;

    MyFrame *frame = new MyFrame("Grid Physics");
    frame->Show(true);
    return true;
}

