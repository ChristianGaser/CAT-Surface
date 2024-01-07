#ifndef  DEF_COLOURS
#define  DEF_COLOURS

/* ----------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1993,1994,1995 David MacDonald,
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */


#include  <volume_io.h>

#define  COMPOSITE_COLOURS( result, front, back ) \
         { \
             int  _rf, _gf, _bf, _af, _rb, _gb, _bb, _ab, _weight; \
             int  _r, _g, _b, _a; \
             _af = get_Colour_a(front); \
             _ab = get_Colour_a(back); \
             if( _af == 255 || _ab == 0 ) \
                 (result) = front; \
             else if( _af == 0 ) \
                 (result) = back; \
             else \
             { \
                 _rf = get_Colour_r(front); \
                 _gf = get_Colour_g(front); \
                 _bf = get_Colour_b(front); \
                 _rb = get_Colour_r(back); \
                 _gb = get_Colour_g(back); \
                 _bb = get_Colour_b(back); \
                 _weight = (255 - _af) * _ab; \
                 _af *= 255; \
                 _r = _af * _rf + _weight * _rb; \
                 _g = _af * _gf + _weight * _gb; \
                 _b = _af * _bf + _weight * _bb; \
                 _a = _af + _weight; \
                 if( _a > 0 ) \
                 { \
                     _r /= _a; \
                     _g /= _a; \
                     _b /= _a; \
                     _a /= 255; \
                 } \
                 (result) = make_rgba_Colour( _r, _g, _b, _a ); \
             } \
         }

typedef  enum  { RGB_SPACE, HSL_SPACE }  Colour_spaces;

#undef  BLACK
#undef  RED
#undef  GREEN
#undef  BLUE
#undef  CYAN
#undef  MAGENTA
#undef  YELLOW
#undef  WHITE

#define ALICE_BLUE make_Colour_0_1( 0.941176, 0.972549, 1.0 )
#define ANTIQUE_WHITE make_Colour_0_1( 0.980392, 0.921569, 0.843137 )
#define AQUAMARINE make_Colour_0_1( 0.498039, 1.0, 0.831373 )
#define AZURE make_Colour_0_1( 0.941176, 1.0, 1.0 )
#define BEIGE make_Colour_0_1( 0.960784, 0.960784, 0.862745 )
#define BISQUE make_Colour_0_1( 1.0, 0.894118, 0.768627 )
#define BLACK make_Colour_0_1( 0.0, 0.0, 0.0 )
#define BLANCHED_ALMOND make_Colour_0_1( 1.0, 0.921569, 0.803922 )
#define BLUE make_Colour_0_1( 0.0, 0.0, 1.0 )
#define BLUE_VIOLET make_Colour_0_1( 0.541176, 0.168627, 0.886275 )
#define BRIGHT_MUSTARD make_Colour_0_1( 1.0, 0.752941, 0.188235 )
#define BROWN make_Colour_0_1( 1.0, 0.25098, 0.25098 )
#define BURLYWOOD make_Colour_0_1( 0.870588, 0.721569, 0.529412 )
#define CADET_BLUE make_Colour_0_1( 0.372549, 0.619608, 0.627451 )
#define CHARTREUSE make_Colour_0_1( 0.498039, 1.0, 0.0 )
#define CHOCOLATE make_Colour_0_1( 0.823529, 0.411765, 0.117647 )
#define CORAL make_Colour_0_1( 1.0, 0.447059, 0.337255 )
#define CORNFLOWER_BLUE make_Colour_0_1( 0.392157, 0.584314, 0.929412 )
#define CORNSILK make_Colour_0_1( 1.0, 0.972549, 0.862745 )
#define CYAN make_Colour_0_1( 0.0, 1.0, 1.0 )
#define DARK_BROWN make_Colour_0_1( 0.313725, 0.156863, 0.117647 )
#define DARK_FIREBRICK make_Colour_0_1( 0.52549, 0.101961, 0.101961 )
#define DARK_GOLDENROD make_Colour_0_1( 1.0, 0.72549, 0.0588235 )
#define DARK_GREEN make_Colour_0_1( 0.0, 0.392157, 0.0 )
#define DARK_KHAKI make_Colour_0_1( 0.741176, 0.717647, 0.419608 )
#define DARK_OLIVE_GREEN make_Colour_0_1( 0.792157, 1.0, 0.439216 )
#define DARK_ORANGE make_Colour_0_1( 1.0, 0.54902, 0.0 )
#define DARK_ORCHID make_Colour_0_1( 0.74902, 0.243137, 1.0 )
#define DARK_RED make_Colour_0_1( 0.501961, 0.0, 0.0 )
#define DARK_SALMON make_Colour_0_1( 0.913725, 0.588235, 0.478431 )
#define DARK_SEA_GREEN make_Colour_0_1( 0.756863, 1.0, 0.756863 )
#define DARK_SLATE_BLUE make_Colour_0_1( 0.282353, 0.239216, 0.545098 )
#define DARK_SLATE_GRAY make_Colour_0_1( 0.592157, 1.0, 1.0 )
#define DARK_SLATE_GREY make_Colour_0_1( 0.184314, 0.309804, 0.309804 )
#define DARK_TURQUOISE make_Colour_0_1( 0.0, 0.807843, 0.819608 )
#define DARK_VIOLET make_Colour_0_1( 0.580392, 0.0, 0.827451 )
#define DEEP_PINK make_Colour_0_1( 1.0, 0.0784314, 0.576471 )
#define DEEP_SKY_BLUE make_Colour_0_1( 0.0, 0.74902, 1.0 )
#define DIM_GRAY make_Colour_0_1( 0.411765, 0.411765, 0.411765 )
#define DIM_GREY make_Colour_0_1( 0.411765, 0.411765, 0.411765 )
#define DODGER_BLUE make_Colour_0_1( 0.117647, 0.564706, 1.0 )
#define FIREBRICK make_Colour_0_1( 0.698039, 0.133333, 0.133333 )
#define FLORAL_WHITE make_Colour_0_1( 1.0, 0.980392, 0.941176 )
#define FOREST_GREEN make_Colour_0_1( 0.133333, 0.545098, 0.133333 )
#define GAINSBORO make_Colour_0_1( 0.862745, 0.862745, 0.862745 )
#define GHOST_WHITE make_Colour_0_1( 0.972549, 0.972549, 1.0 )
#define GOLDENROD make_Colour_0_1( 1.0, 0.756863, 0.145098 )
#define GOLD make_Colour_0_1( 1.0, 0.843137, 0.0 )
#define GRAY make_Colour_0_1( 0.752941, 0.752941, 0.752941 )
#define GREEN make_Colour_0_1( 0.0, 1.0, 0.0 )
#define GREEN_YELLOW make_Colour_0_1( 0.678431, 1.0, 0.184314 )
#define GREY make_Colour_0_1( 0.752941, 0.752941, 0.752941 )
#define HONEYDEW make_Colour_0_1( 0.941176, 1.0, 0.941176 )
#define HOT_PINK make_Colour_0_1( 1.0, 0.431373, 0.705882 )
#define INDIAN_RED make_Colour_0_1( 1.0, 0.415686, 0.415686 )
#define IVORY make_Colour_0_1( 1.0, 1.0, 0.941176 )
#define KHAKI make_Colour_0_1( 1.0, 0.964706, 0.560784 )
#define LAVENDER_BLUSH make_Colour_0_1( 1.0, 0.941176, 0.960784 )
#define LAVENDER make_Colour_0_1( 0.901961, 0.901961, 0.980392 )
#define LAWN_GREEN make_Colour_0_1( 0.486275, 0.988235, 0.0 )
#define LEMON_CHIFFON make_Colour_0_1( 1.0, 0.980392, 0.803922 )
#define LIGHT_BLUE make_Colour_0_1( 0.74902, 0.937255, 1.0 )
#define LIGHT_CORAL make_Colour_0_1( 0.941176, 0.501961, 0.501961 )
#define LIGHT_CYAN make_Colour_0_1( 0.878431, 1.0, 1.0 )
#define LIGHT_GOLDENROD make_Colour_0_1( 1.0, 0.92549, 0.545098 )
#define LIGHT_GOLDENROD_YELLOW make_Colour_0_1( 0.980392, 0.980392, 0.823529 )
#define LIGHT_GRAY make_Colour_0_1( 0.827451, 0.827451, 0.827451 )
#define LIGHT_GREY make_Colour_0_1( 0.827451, 0.827451, 0.827451 )
#define LIGHT_PINK make_Colour_0_1( 1.0, 0.713725, 0.756863 )
#define LIGHT_SALMON make_Colour_0_1( 1.0, 0.627451, 0.478431 )
#define LIGHT_SEA_GREEN make_Colour_0_1( 0.12549, 0.698039, 0.666667 )
#define LIGHT_SKY_BLUE make_Colour_0_1( 0.529412, 0.807843, 0.980392 )
#define LIGHT_SLATE_BLUE make_Colour_0_1( 0.517647, 0.439216, 1.0 )
#define LIGHT_SLATE_GRAY make_Colour_0_1( 0.466667, 0.533333, 0.6 )
#define LIGHT_SLATE_GREY make_Colour_0_1( 0.466667, 0.533333, 0.6 )
#define LIGHT_STEEL_BLUE make_Colour_0_1( 0.792157, 0.882353, 1.0 )
#define LIGHT_YELLOW make_Colour_0_1( 1.0, 1.0, 0.878431 )
#define LIME_GREEN make_Colour_0_1( 0.196078, 0.803922, 0.196078 )
#define LINEN make_Colour_0_1( 0.980392, 0.941176, 0.901961 )
#define MAGENTA make_Colour_0_1( 1.0, 0.0, 1.0 )
#define MAROON make_Colour_0_1( 1.0, 0.203922, 0.701961 )
#define MEDIUM_AQUAMARINE make_Colour_0_1( 0.4, 0.803922, 0.666667 )
#define MEDIUM_BLUE make_Colour_0_1( 0.0, 0.0, 0.803922 )
#define MEDIUM_ORCHID make_Colour_0_1( 0.729412, 0.333333, 0.827451 )
#define MEDIUM_PURPLE make_Colour_0_1( 0.670588, 0.509804, 1.0 )
#define MEDIUM_SEA_GREEN make_Colour_0_1( 0.235294, 0.701961, 0.443137 )
#define MEDIUM_SLATE_BLUE make_Colour_0_1( 0.482353, 0.407843, 0.933333 )
#define MEDIUM_SPRING_GREEN make_Colour_0_1( 0.0, 0.980392, 0.603922 )
#define MEDIUM_TURQUOISE make_Colour_0_1( 0.282353, 0.819608, 0.8 )
#define MEDIUM_VIOLET_RED make_Colour_0_1( 0.780392, 0.0823529, 0.521569 )
#define MIAMI_PINK make_Colour_0_1( 1.0, 0.501961, 0.564706 )
#define MIAMI_TURQUOISE make_Colour_0_1( 0.313725, 1.0, 1.0 )
#define MIDNIGHT_BLUE make_Colour_0_1( 0.0980392, 0.0980392, 0.439216 )
#define MINT_CREAM make_Colour_0_1( 0.960784, 1.0, 0.980392 )
#define MISTY_ROSE make_Colour_0_1( 1.0, 0.894118, 0.882353 )
#define MOCCASIN make_Colour_0_1( 1.0, 0.894118, 0.709804 )
#define NAVAJO_WHITE make_Colour_0_1( 1.0, 0.870588, 0.678431 )
#define NAVY_BLUE make_Colour_0_1( 0.0, 0.0, 0.501961 )
#define NAVY make_Colour_0_1( 0.0, 0.0, 0.501961 )
#define OLD_LACE make_Colour_0_1( 0.992157, 0.960784, 0.901961 )
#define OLIVE_DRAB make_Colour_0_1( 0.419608, 0.556863, 0.137255 )
#define ORANGE make_Colour_0_1( 1.0, 0.647059, 0.0 )
#define ORANGE_RED make_Colour_0_1( 1.0, 0.270588, 0.0 )
#define ORCHID make_Colour_0_1( 0.854902, 0.439216, 0.839216 )
#define PALE_GOLDENROD make_Colour_0_1( 0.933333, 0.909804, 0.666667 )
#define PALE_GREEN make_Colour_0_1( 0.603922, 1.0, 0.603922 )
#define PALE_TURQUOISE make_Colour_0_1( 0.733333, 1.0, 1.0 )
#define PALE_VIOLET_RED make_Colour_0_1( 1.0, 0.509804, 0.670588 )
#define PAPAYA_WHIP make_Colour_0_1( 1.0, 0.937255, 0.835294 )
#define PEACH_PUFF make_Colour_0_1( 1.0, 0.854902, 0.72549 )
#define PERU make_Colour_0_1( 0.803922, 0.521569, 0.247059 )
#define PINK make_Colour_0_1( 1.0, 0.752941, 0.796078 )
#define PLUM make_Colour_0_1( 1.0, 0.733333, 1.0 )
#define POWDER_BLUE make_Colour_0_1( 0.690196, 0.878431, 0.901961 )
#define PURPLE make_Colour_0_1( 0.627451, 0.12549, 0.941176 )
#define RED make_Colour_0_1( 1.0, 0.0, 0.0 )
#define ROSY_BROWN make_Colour_0_1( 1.0, 0.756863, 0.756863 )
#define ROYAL_BLUE make_Colour_0_1( 0.282353, 0.462745, 1.0 )
#define SADDLE_BROWN make_Colour_0_1( 0.545098, 0.270588, 0.0745098 )
#define SALMON make_Colour_0_1( 1.0, 0.54902, 0.411765 )
#define SANDY_BROWN make_Colour_0_1( 0.956863, 0.643137, 0.376471 )
#define SEASHELL make_Colour_0_1( 1.0, 0.960784, 0.933333 )
#define SEA_GREEN make_Colour_0_1( 0.329412, 1.0, 0.623529 )
#define SIENNA make_Colour_0_1( 1.0, 0.509804, 0.278431 )
#define SKY_BLUE make_Colour_0_1( 0.529412, 0.807843, 1.0 )
#define SLATE_BLUE make_Colour_0_1( 0.513725, 0.435294, 1.0 )
#define SLATE_GRAY make_Colour_0_1( 0.776471, 0.886275, 1.0 )
#define SLATE_GREY make_Colour_0_1( 0.439216, 0.501961, 0.564706 )
#define SNOW make_Colour_0_1( 1.0, 0.980392, 0.980392 )
#define SPRING_GREEN make_Colour_0_1( 0.0, 1.0, 0.498039 )
#define STEEL_BLUE make_Colour_0_1( 0.388235, 0.721569, 1.0 )
#define TAN make_Colour_0_1( 1.0, 0.647059, 0.309804 )
#define THISTLE make_Colour_0_1( 1.0, 0.882353, 1.0 )
#define TOMATO make_Colour_0_1( 1.0, 0.388235, 0.278431 )
#define TURQUOISE make_Colour_0_1( 0.25098, 0.878431, 0.815686 )
#define VIOLET make_Colour_0_1( 0.933333, 0.509804, 0.933333 )
#define VIOLET_RED make_Colour_0_1( 1.0, 0.243137, 0.588235 )
#define WHEAT make_Colour_0_1( 1.0, 0.905882, 0.729412 )
#define WHITE make_Colour_0_1( 1.0, 1.0, 1.0 )
#define WHITE_SMOKE make_Colour_0_1( 0.960784, 0.960784, 0.960784 )
#define YELLOW make_Colour_0_1( 1.0, 1.0, 0.0 )
#define YELLOW_GREEN make_Colour_0_1( 0.603922, 0.803922, 0.196078 )

#endif
