#ifndef __PANGOVERWRITE__
#define __PANGOVERWRITE__

#include <pangolin/display/display_internal.h>
#include <pangolin/gl/glfont.h>

//Overwrite some of pangolin's internal functions.
extern "C" const unsigned char AnonymousPro_ttf[];

namespace pangolin
{

    template <typename T>
    void GuiVarChanged(Var<T> &var)
    {
        VarState::I().FlagVarChanged();
        var.Meta().gui_changed = true;

        for (std::vector<GuiVarChangedCallback>::iterator igvc = VarState::I().gui_var_changed_callbacks.begin(); igvc != VarState::I().gui_var_changed_callbacks.end(); ++igvc)
        {
            if (StartsWith(var.Meta().full_name, igvc->filter))
            {
                igvc->fn(igvc->data, var.Meta().full_name, var.Ref());
            }
        }
    }

    std::mutex new_display_mutex;

    void glRect(Viewport v)
    {
        GLfloat vs[] = {(float)v.l, (float)v.b,
                        (float)v.l, (float)v.t(),
                        (float)v.r(), (float)v.t(),
                        (float)v.r(), (float)v.b};

        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(2, GL_FLOAT, 0, vs);
        glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
        glDisableClientState(GL_VERTEX_ARRAY);
    }

    struct PANGOLIN_EXPORT HandlerResize : Handler
    {

        void Mouse(View &In_, MouseButton button, int x, int y, bool pressed, int button_state)
        {
            float TopBound = In_.top.p;

            if (button == MouseWheelDown) //MouseButtonRight
                TopBound = std::max(0.1, TopBound - 0.01);

            if (button == MouseWheelUp) //MouseButtonLeft
                TopBound = std::min(0.5, TopBound + 0.01);

            In_.SetBounds(In_.bottom, TopBound, In_.left, In_.right);
        }

        void MouseMotion(View &, int x, int y, int button_state)
        {
        }
    };

} // namespace pangolin

#endif
