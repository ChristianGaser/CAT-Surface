
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

#include "bicpl_internal.h"

#ifndef lint
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Objects/colours.c,v 1.11 2005/08/17 22:28:26 bert Exp $";
#endif

static Colour ANTIQUE_WHITE_COL;
static Colour ALICE_BLUE_COL;
static Colour AQUAMARINE_COL;
static Colour AZURE_COL;
static Colour BEIGE_COL;
static Colour BISQUE_COL;
static Colour BLACK_COL;
static Colour BLANCHED_ALMOND_COL;
static Colour BLUE_COL;
static Colour BLUE_VIOLET_COL;
static Colour BRIGHT_MUSTARD_COL;
static Colour BROWN_COL;
static Colour BURLYWOOD_COL;
static Colour CADET_BLUE_COL;
static Colour CHARTREUSE_COL;
static Colour CHOCOLATE_COL;
static Colour CORAL_COL;
static Colour CORNFLOWER_BLUE_COL;
static Colour CORNSILK_COL;
static Colour CYAN_COL;
static Colour DARK_BROWN_COL;
static Colour DARK_FIREBRICK_COL;
static Colour DARK_GOLDENROD_COL;
static Colour DARK_GREEN_COL;
static Colour DARK_KHAKI_COL;
static Colour DARK_OLIVE_GREEN_COL;
static Colour DARK_ORANGE_COL;
static Colour DARK_ORCHID_COL;
static Colour DARK_RED_COL;
static Colour DARK_SALMON_COL;
static Colour DARK_SEA_GREEN_COL;
static Colour DARK_SLATE_BLUE_COL;
static Colour DARK_SLATE_GRAY_COL;
static Colour DARK_SLATE_GREY_COL;
static Colour DARK_TURQUOISE_COL;
static Colour DARK_VIOLET_COL;
static Colour DEEP_PINK_COL;
static Colour DEEP_SKY_BLUE_COL;
static Colour DIM_GRAY_COL;
static Colour DIM_GREY_COL;
static Colour DODGER_BLUE_COL;
static Colour FIREBRICK_COL;
static Colour FLORAL_WHITE_COL;
static Colour FOREST_GREEN_COL;
static Colour GAINSBORO_COL;
static Colour GHOST_WHITE_COL;
static Colour GOLDENROD_COL;
static Colour GOLD_COL;
static Colour GRAY_COL;
static Colour GREEN_COL;
static Colour GREEN_YELLOW_COL;
static Colour GREY_COL;
static Colour HONEYDEW_COL;
static Colour HOT_PINK_COL;
static Colour INDIAN_RED_COL;
static Colour IVORY_COL;
static Colour KHAKI_COL;
static Colour LAVENDER_BLUSH_COL;
static Colour LAVENDER_COL;
static Colour LAWN_GREEN_COL;
static Colour LEMON_CHIFFON_COL;
static Colour LIGHT_BLUE_COL;
static Colour LIGHT_CORAL_COL;
static Colour LIGHT_CYAN_COL;
static Colour LIGHT_GOLDENROD_COL;
static Colour LIGHT_GOLDENROD_YELLOW_COL;
static Colour LIGHT_GRAY_COL;
static Colour LIGHT_GREY_COL;
static Colour LIGHT_PINK_COL;
static Colour LIGHT_SALMON_COL;
static Colour LIGHT_SEA_GREEN_COL;
static Colour LIGHT_SKY_BLUE_COL;
static Colour LIGHT_SLATE_BLUE_COL;
static Colour LIGHT_SLATE_GRAY_COL;
static Colour LIGHT_SLATE_GREY_COL;
static Colour LIGHT_STEEL_BLUE_COL;
static Colour LIGHT_YELLOW_COL;
static Colour LIME_GREEN_COL;
static Colour LINEN_COL;
static Colour MAGENTA_COL;
static Colour MAROON_COL;
static Colour MEDIUM_AQUAMARINE_COL;
static Colour MEDIUM_BLUE_COL;
static Colour MEDIUM_ORCHID_COL;
static Colour MEDIUM_PURPLE_COL;
static Colour MEDIUM_SEA_GREEN_COL;
static Colour MEDIUM_SLATE_BLUE_COL;
static Colour MEDIUM_SPRING_GREEN_COL;
static Colour MEDIUM_TURQUOISE_COL;
static Colour MEDIUM_VIOLET_RED_COL;
static Colour MIAMI_PINK_COL;
static Colour MIAMI_TURQUOISE_COL;
static Colour MIDNIGHT_BLUE_COL;
static Colour MINT_CREAM_COL;
static Colour MISTY_ROSE_COL;
static Colour MOCCASIN_COL;
static Colour NAVAJO_WHITE_COL;
static Colour NAVY_BLUE_COL;
static Colour NAVY_COL;
static Colour OLD_LACE_COL;
static Colour OLIVE_DRAB_COL;
static Colour ORANGE_COL;
static Colour ORANGE_RED_COL;
static Colour ORCHID_COL;
static Colour PALE_GOLDENROD_COL;
static Colour PALE_GREEN_COL;
static Colour PALE_TURQUOISE_COL;
static Colour PALE_VIOLET_RED_COL;
static Colour PAPAYA_WHIP_COL;
static Colour PEACH_PUFF_COL;
static Colour PERU_COL;
static Colour PINK_COL;
static Colour PLUM_COL;
static Colour POWDER_BLUE_COL;
static Colour PURPLE_COL;
static Colour RED_COL;
static Colour ROSY_BROWN_COL;
static Colour ROYAL_BLUE_COL;
static Colour SADDLE_BROWN_COL;
static Colour SALMON_COL;
static Colour SANDY_BROWN_COL;
static Colour SEASHELL_COL;
static Colour SEA_GREEN_COL;
static Colour SIENNA_COL;
static Colour SKY_BLUE_COL;
static Colour SLATE_BLUE_COL;
static Colour SLATE_GRAY_COL;
static Colour SLATE_GREY_COL;
static Colour SNOW_COL;
static Colour SPRING_GREEN_COL;
static Colour STEEL_BLUE_COL;
static Colour TAN_COL;
static Colour THISTLE_COL;
static Colour TOMATO_COL;
static Colour TURQUOISE_COL;
static Colour VIOLET_COL;
static Colour VIOLET_RED_COL;
static Colour WHEAT_COL;
static Colour WHITE_COL;
static Colour WHITE_SMOKE_COL;
static Colour YELLOW_COL;
static Colour YELLOW_GREEN_COL;
static Colour TRANSPARENT_COL;

typedef  struct
{
    STRING   name;
    Colour   *colour;
} colours_struct;

static  colours_struct  colour_lookup[] =
{
    { "TRANSPARENT", &TRANSPARENT_COL },
    { "ALICE_BLUE", &ALICE_BLUE_COL },
    { "ANTIQUE_WHITE", &ANTIQUE_WHITE_COL },
    { "AQUAMARINE", &AQUAMARINE_COL },
    { "AZURE", &AZURE_COL },
    { "BEIGE", &BEIGE_COL },
    { "BISQUE", &BISQUE_COL },
    { "BLACK", &BLACK_COL },
    { "BLANCHED_ALMOND", &BLANCHED_ALMOND_COL },
    { "BLUE", &BLUE_COL },
    { "BLUE_VIOLET", &BLUE_VIOLET_COL },
    { "BRIGHT_MUSTARD", &BRIGHT_MUSTARD_COL },
    { "BROWN", &BROWN_COL },
    { "BURLYWOOD", &BURLYWOOD_COL },
    { "CADET_BLUE", &CADET_BLUE_COL },
    { "CHARTREUSE", &CHARTREUSE_COL },
    { "CHOCOLATE", &CHOCOLATE_COL },
    { "CORAL", &CORAL_COL },
    { "CORNFLOWER_BLUE", &CORNFLOWER_BLUE_COL },
    { "CORNSILK", &CORNSILK_COL },
    { "CYAN", &CYAN_COL },
    { "DARK_BROWN", &DARK_BROWN_COL },
    { "DARK_FIREBRICK", &DARK_FIREBRICK_COL },
    { "DARK_GOLDENROD", &DARK_GOLDENROD_COL },
    { "DARK_GREEN", &DARK_GREEN_COL },
    { "DARK_KHAKI", &DARK_KHAKI_COL },
    { "DARK_OLIVE_GREEN", &DARK_OLIVE_GREEN_COL },
    { "DARK_ORANGE", &DARK_ORANGE_COL },
    { "DARK_ORCHID", &DARK_ORCHID_COL },
    { "DARK_RED", &DARK_RED_COL },
    { "DARK_SALMON", &DARK_SALMON_COL },
    { "DARK_SEA_GREEN", &DARK_SEA_GREEN_COL },
    { "DARK_SLATE_BLUE", &DARK_SLATE_BLUE_COL },
    { "DARK_SLATE_GRAY", &DARK_SLATE_GRAY_COL },
    { "DARK_SLATE_GREY", &DARK_SLATE_GREY_COL },
    { "DARK_TURQUOISE", &DARK_TURQUOISE_COL },
    { "DARK_VIOLET", &DARK_VIOLET_COL },
    { "DEEP_PINK", &DEEP_PINK_COL },
    { "DEEP_SKY_BLUE", &DEEP_SKY_BLUE_COL },
    { "DIM_GRAY", &DIM_GRAY_COL },
    { "DIM_GREY", &DIM_GREY_COL },
    { "DODGER_BLUE", &DODGER_BLUE_COL },
    { "FIREBRICK", &FIREBRICK_COL },
    { "FLORAL_WHITE", &FLORAL_WHITE_COL },
    { "FOREST_GREEN", &FOREST_GREEN_COL },
    { "GAINSBORO", &GAINSBORO_COL },
    { "GHOST_WHITE", &GHOST_WHITE_COL },
    { "GOLD", &GOLD_COL },
    { "GOLDENROD", &GOLDENROD_COL },
    { "GRAY", &GRAY_COL },
    { "GREEN", &GREEN_COL },
    { "GREEN_YELLOW", &GREEN_YELLOW_COL },
    { "GREY", &GREY_COL },
    { "HONEYDEW", &HONEYDEW_COL },
    { "HOT_PINK", &HOT_PINK_COL },
    { "INDIAN_RED", &INDIAN_RED_COL },
    { "IVORY", &IVORY_COL },
    { "KHAKI", &KHAKI_COL },
    { "LAVENDER", &LAVENDER_COL },
    { "LAVENDER_BLUSH", &LAVENDER_BLUSH_COL },
    { "LAWN_GREEN", &LAWN_GREEN_COL },
    { "LEMON_CHIFFON", &LEMON_CHIFFON_COL },
    { "LIGHT_BLUE", &LIGHT_BLUE_COL },
    { "LIGHT_CORAL", &LIGHT_CORAL_COL },
    { "LIGHT_CYAN", &LIGHT_CYAN_COL },
    { "LIGHT_GOLDENROD", &LIGHT_GOLDENROD_COL },
    { "LIGHT_GOLDENROD_YELLOW", &LIGHT_GOLDENROD_YELLOW_COL },
    { "LIGHT_GRAY", &LIGHT_GRAY_COL },
    { "LIGHT_GREY", &LIGHT_GREY_COL },
    { "LIGHT_PINK", &LIGHT_PINK_COL },
    { "LIGHT_SALMON", &LIGHT_SALMON_COL },
    { "LIGHT_SEA_GREEN", &LIGHT_SEA_GREEN_COL },
    { "LIGHT_SKY_BLUE", &LIGHT_SKY_BLUE_COL },
    { "LIGHT_SLATE_BLUE", &LIGHT_SLATE_BLUE_COL },
    { "LIGHT_SLATE_GRAY", &LIGHT_SLATE_GRAY_COL },
    { "LIGHT_SLATE_GREY", &LIGHT_SLATE_GREY_COL },
    { "LIGHT_STEEL_BLUE", &LIGHT_STEEL_BLUE_COL },
    { "LIGHT_YELLOW", &LIGHT_YELLOW_COL },
    { "LIME_GREEN", &LIME_GREEN_COL },
    { "LINEN", &LINEN_COL },
    { "MAGENTA", &MAGENTA_COL },
    { "MAROON", &MAROON_COL },
    { "MEDIUM_AQUAMARINE", &MEDIUM_AQUAMARINE_COL },
    { "MEDIUM_BLUE", &MEDIUM_BLUE_COL },
    { "MEDIUM_ORCHID", &MEDIUM_ORCHID_COL },
    { "MEDIUM_PURPLE", &MEDIUM_PURPLE_COL },
    { "MEDIUM_SEA_GREEN", &MEDIUM_SEA_GREEN_COL },
    { "MEDIUM_SLATE_BLUE", &MEDIUM_SLATE_BLUE_COL },
    { "MEDIUM_SPRING_GREEN", &MEDIUM_SPRING_GREEN_COL },
    { "MEDIUM_TURQUOISE", &MEDIUM_TURQUOISE_COL },
    { "MEDIUM_VIOLET_RED", &MEDIUM_VIOLET_RED_COL },
    { "MIAMI_PINK", &MIAMI_PINK_COL },
    { "MIAMI_TURQUOISE", &MIAMI_TURQUOISE_COL },
    { "MIDNIGHT_BLUE", &MIDNIGHT_BLUE_COL },
    { "MINT_CREAM", &MINT_CREAM_COL },
    { "MISTY_ROSE", &MISTY_ROSE_COL },
    { "MOCCASIN", &MOCCASIN_COL },
    { "NAVAJO_WHITE", &NAVAJO_WHITE_COL },
    { "NAVY", &NAVY_COL },
    { "NAVY_BLUE", &NAVY_BLUE_COL },
    { "OLD_LACE", &OLD_LACE_COL },
    { "OLIVE_DRAB", &OLIVE_DRAB_COL },
    { "ORANGE", &ORANGE_COL },
    { "ORANGE_RED", &ORANGE_RED_COL },
    { "ORCHID", &ORCHID_COL },
    { "PALE_GOLDENROD", &PALE_GOLDENROD_COL },
    { "PALE_GREEN", &PALE_GREEN_COL },
    { "PALE_TURQUOISE", &PALE_TURQUOISE_COL },
    { "PALE_VIOLET_RED", &PALE_VIOLET_RED_COL },
    { "PAPAYA_WHIP", &PAPAYA_WHIP_COL },
    { "PEACH_PUFF", &PEACH_PUFF_COL },
    { "PERU", &PERU_COL },
    { "PINK", &PINK_COL },
    { "PLUM", &PLUM_COL },
    { "POWDER_BLUE", &POWDER_BLUE_COL },
    { "PURPLE", &PURPLE_COL },
    { "RED", &RED_COL },
    { "ROSY_BROWN", &ROSY_BROWN_COL },
    { "ROYAL_BLUE", &ROYAL_BLUE_COL },
    { "SADDLE_BROWN", &SADDLE_BROWN_COL },
    { "SALMON", &SALMON_COL },
    { "SANDY_BROWN", &SANDY_BROWN_COL },
    { "SEASHELL", &SEASHELL_COL },
    { "SEA_GREEN", &SEA_GREEN_COL },
    { "SIENNA", &SIENNA_COL },
    { "SKY_BLUE", &SKY_BLUE_COL },
    { "SLATE_BLUE", &SLATE_BLUE_COL },
    { "SLATE_GRAY", &SLATE_GRAY_COL },
    { "SLATE_GREY", &SLATE_GREY_COL },
    { "SNOW", &SNOW_COL },
    { "SPRING_GREEN", &SPRING_GREEN_COL },
    { "STEEL_BLUE", &STEEL_BLUE_COL },
    { "TAN", &TAN_COL },
    { "THISTLE", &THISTLE_COL },
    { "TOMATO", &TOMATO_COL },
    { "TURQUOISE", &TURQUOISE_COL },
    { "VIOLET", &VIOLET_COL },
    { "VIOLET_RED", &VIOLET_RED_COL },
    { "WHEAT", &WHEAT_COL },
    { "WHITE", &WHITE_COL },
    { "WHITE_SMOKE", &WHITE_SMOKE_COL },
    { "YELLOW", &YELLOW_COL },
    { "YELLOW_GREEN", &YELLOW_GREEN_COL }
};

static   BOOLEAN  strings_equivalent( STRING, STRING );

/* ----------------------------- MNI Header -----------------------------------
@NAME       : check_initialize_colours
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: The first time it is called, this function creates the global
              colour variables.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  check_initialize_colours( void )
{
    static  BOOLEAN  first = TRUE;

    if( first )
    {
        first = FALSE;

        ANTIQUE_WHITE_COL = make_Colour_0_1( 0.980392, 0.921569, 0.843137 );
        AQUAMARINE_COL = make_Colour_0_1( 0.498039, 1.0, 0.831373 );
        AZURE_COL = make_Colour_0_1( 0.941176, 1.0, 1.0 );
        BEIGE_COL = make_Colour_0_1( 0.960784, 0.960784, 0.862745 );
        BISQUE_COL = make_Colour_0_1( 1.0, 0.894118, 0.768627 );
        BLACK_COL = make_Colour_0_1( 0.0, 0.0, 0.0 );
        BLANCHED_ALMOND_COL = make_Colour_0_1( 1.0, 0.921569, 0.803922 );
        BLUE_COL = make_Colour_0_1( 0.0, 0.0, 1.0 );
        BLUE_VIOLET_COL = make_Colour_0_1( 0.541176, 0.168627, 0.886275 );
        BRIGHT_MUSTARD_COL = make_Colour_0_1( 1.0, 0.752941, 0.188235 );
        BROWN_COL = make_Colour_0_1( 1.0, 0.25098, 0.25098 );
        BURLYWOOD_COL = make_Colour_0_1( 0.870588, 0.721569, 0.529412 );
        CADET_BLUE_COL = make_Colour_0_1( 0.372549, 0.619608, 0.627451 );
        CHARTREUSE_COL = make_Colour_0_1( 0.498039, 1.0, 0.0 );
        CHOCOLATE_COL = make_Colour_0_1( 0.823529, 0.411765, 0.117647 );
        CORAL_COL = make_Colour_0_1( 1.0, 0.447059, 0.337255 );
        CORNFLOWER_BLUE_COL = make_Colour_0_1( 0.392157, 0.584314, 0.929412 );
        CORNSILK_COL = make_Colour_0_1( 1.0, 0.972549, 0.862745 );
        CYAN_COL = make_Colour_0_1( 0.0, 1.0, 1.0 );
        DARK_BROWN_COL = make_Colour_0_1( 0.313725, 0.156863, 0.117647 );
        DARK_FIREBRICK_COL = make_Colour_0_1( 0.52549, 0.101961, 0.101961 );
        DARK_GOLDENROD_COL = make_Colour_0_1( 1.0, 0.72549, 0.0588235 );
        DARK_GREEN_COL = make_Colour_0_1( 0.0, 0.392157, 0.0 );
        DARK_KHAKI_COL = make_Colour_0_1( 0.741176, 0.717647, 0.419608 );
        DARK_OLIVE_GREEN_COL = make_Colour_0_1( 0.792157, 1.0, 0.439216 );
        DARK_ORANGE_COL = make_Colour_0_1( 1.0, 0.54902, 0.0 );
        DARK_ORCHID_COL = make_Colour_0_1( 0.74902, 0.243137, 1.0 );
        DARK_RED_COL = make_Colour_0_1( 0.501961, 0.0, 0.0 );
        DARK_SALMON_COL = make_Colour_0_1( 0.913725, 0.588235, 0.478431 );
        DARK_SEA_GREEN_COL = make_Colour_0_1( 0.756863, 1.0, 0.756863 );
        DARK_SLATE_BLUE_COL = make_Colour_0_1( 0.282353, 0.239216, 0.545098 );
        DARK_SLATE_GRAY_COL = make_Colour_0_1( 0.592157, 1.0, 1.0 );
        DARK_SLATE_GREY_COL = make_Colour_0_1( 0.184314, 0.309804, 0.309804 );
        DARK_TURQUOISE_COL = make_Colour_0_1( 0.0, 0.807843, 0.819608 );
        DARK_VIOLET_COL = make_Colour_0_1( 0.580392, 0.0, 0.827451 );
        DEEP_PINK_COL = make_Colour_0_1( 1.0, 0.0784314, 0.576471 );
        DEEP_SKY_BLUE_COL = make_Colour_0_1( 0.0, 0.74902, 1.0 );
        DIM_GRAY_COL = make_Colour_0_1( 0.411765, 0.411765, 0.411765 );
        DIM_GREY_COL = make_Colour_0_1( 0.411765, 0.411765, 0.411765 );
        DODGER_BLUE_COL = make_Colour_0_1( 0.117647, 0.564706, 1.0 );
        FIREBRICK_COL = make_Colour_0_1( 0.698039, 0.133333, 0.133333 );
        FLORAL_WHITE_COL = make_Colour_0_1( 1.0, 0.980392, 0.941176 );
        FOREST_GREEN_COL = make_Colour_0_1( 0.133333, 0.545098, 0.133333 );
        GAINSBORO_COL = make_Colour_0_1( 0.862745, 0.862745, 0.862745 );
        GHOST_WHITE_COL = make_Colour_0_1( 0.972549, 0.972549, 1.0 );
        GOLDENROD_COL = make_Colour_0_1( 1.0, 0.756863, 0.145098 );
        GOLD_COL = make_Colour_0_1( 1.0, 0.843137, 0.0 );
        GRAY_COL = make_Colour_0_1( 0.752941, 0.752941, 0.752941 );
        GREEN_COL = make_Colour_0_1( 0.0, 1.0, 0.0 );
        GREEN_YELLOW_COL = make_Colour_0_1( 0.678431, 1.0, 0.184314 );
        GREY_COL = make_Colour_0_1( 0.752941, 0.752941, 0.752941 );
        HONEYDEW_COL = make_Colour_0_1( 0.941176, 1.0, 0.941176 );
        HOT_PINK_COL = make_Colour_0_1( 1.0, 0.431373, 0.705882 );
        INDIAN_RED_COL = make_Colour_0_1( 1.0, 0.415686, 0.415686 );
        IVORY_COL = make_Colour_0_1( 1.0, 1.0, 0.941176 );
        KHAKI_COL = make_Colour_0_1( 1.0, 0.964706, 0.560784 );
        LAVENDER_BLUSH_COL = make_Colour_0_1( 1.0, 0.941176, 0.960784 );
        LAVENDER_COL = make_Colour_0_1( 0.901961, 0.901961, 0.980392 );
        LAWN_GREEN_COL = make_Colour_0_1( 0.486275, 0.988235, 0.0 );
        LEMON_CHIFFON_COL = make_Colour_0_1( 1.0, 0.980392, 0.803922 );
        LIGHT_BLUE_COL = make_Colour_0_1( 0.74902, 0.937255, 1.0 );
        LIGHT_CORAL_COL = make_Colour_0_1( 0.941176, 0.501961, 0.501961 );
        LIGHT_CYAN_COL = make_Colour_0_1( 0.878431, 1.0, 1.0 );
        LIGHT_GOLDENROD_COL = make_Colour_0_1( 1.0, 0.92549, 0.545098 );
        LIGHT_GOLDENROD_YELLOW_COL = make_Colour_0_1( 0.980392, 0.980392, 0.823529 );
        LIGHT_GRAY_COL = make_Colour_0_1( 0.827451, 0.827451, 0.827451 );
        LIGHT_GREY_COL = make_Colour_0_1( 0.827451, 0.827451, 0.827451 );
        LIGHT_PINK_COL = make_Colour_0_1( 1.0, 0.713725, 0.756863 );
        LIGHT_SALMON_COL = make_Colour_0_1( 1.0, 0.627451, 0.478431 );
        LIGHT_SEA_GREEN_COL = make_Colour_0_1( 0.12549, 0.698039, 0.666667 );
        LIGHT_SKY_BLUE_COL = make_Colour_0_1( 0.529412, 0.807843, 0.980392 );
        LIGHT_SLATE_BLUE_COL = make_Colour_0_1( 0.517647, 0.439216, 1.0 );
        LIGHT_SLATE_GRAY_COL = make_Colour_0_1( 0.466667, 0.533333, 0.6 );
        LIGHT_SLATE_GREY_COL = make_Colour_0_1( 0.466667, 0.533333, 0.6 );
        LIGHT_STEEL_BLUE_COL = make_Colour_0_1( 0.792157, 0.882353, 1.0 );
        LIGHT_YELLOW_COL = make_Colour_0_1( 1.0, 1.0, 0.878431 );
        LIME_GREEN_COL = make_Colour_0_1( 0.196078, 0.803922, 0.196078 );
        LINEN_COL = make_Colour_0_1( 0.980392, 0.941176, 0.901961 );
        MAGENTA_COL = make_Colour_0_1( 1.0, 0.0, 1.0 );
        MAROON_COL = make_Colour_0_1( 1.0, 0.203922, 0.701961 );
        MEDIUM_AQUAMARINE_COL = make_Colour_0_1( 0.4, 0.803922, 0.666667 );
        MEDIUM_BLUE_COL = make_Colour_0_1( 0.0, 0.0, 0.803922 );
        MEDIUM_ORCHID_COL = make_Colour_0_1( 0.729412, 0.333333, 0.827451 );
        MEDIUM_PURPLE_COL = make_Colour_0_1( 0.670588, 0.509804, 1.0 );
        MEDIUM_SEA_GREEN_COL = make_Colour_0_1( 0.235294, 0.701961, 0.443137 );
        MEDIUM_SLATE_BLUE_COL = make_Colour_0_1( 0.482353, 0.407843, 0.933333 );
        MEDIUM_SPRING_GREEN_COL = make_Colour_0_1( 0.0, 0.980392, 0.603922 );
        MEDIUM_TURQUOISE_COL = make_Colour_0_1( 0.282353, 0.819608, 0.8 );
        MEDIUM_VIOLET_RED_COL = make_Colour_0_1( 0.780392, 0.0823529, 0.521569 );
        MIAMI_PINK_COL = make_Colour_0_1( 1.0, 0.501961, 0.564706 );
        MIAMI_TURQUOISE_COL = make_Colour_0_1( 0.313725, 1.0, 1.0 );
        MIDNIGHT_BLUE_COL = make_Colour_0_1( 0.0980392, 0.0980392, 0.439216 );
        MINT_CREAM_COL = make_Colour_0_1( 0.960784, 1.0, 0.980392 );
        MISTY_ROSE_COL = make_Colour_0_1( 1.0, 0.894118, 0.882353 );
        MOCCASIN_COL = make_Colour_0_1( 1.0, 0.894118, 0.709804 );
        NAVAJO_WHITE_COL = make_Colour_0_1( 1.0, 0.870588, 0.678431 );
        NAVY_BLUE_COL = make_Colour_0_1( 0.0, 0.0, 0.501961 );
        NAVY_COL = make_Colour_0_1( 0.0, 0.0, 0.501961 );
        OLD_LACE_COL = make_Colour_0_1( 0.992157, 0.960784, 0.901961 );
        OLIVE_DRAB_COL = make_Colour_0_1( 0.419608, 0.556863, 0.137255 );
        ORANGE_COL = make_Colour_0_1( 1.0, 0.647059, 0.0 );
        ORANGE_RED_COL = make_Colour_0_1( 1.0, 0.270588, 0.0 );
        ORCHID_COL = make_Colour_0_1( 0.854902, 0.439216, 0.839216 );
        PALE_GOLDENROD_COL = make_Colour_0_1( 0.933333, 0.909804, 0.666667 );
        PALE_GREEN_COL = make_Colour_0_1( 0.603922, 1.0, 0.603922 );
        PALE_TURQUOISE_COL = make_Colour_0_1( 0.733333, 1.0, 1.0 );
        PALE_VIOLET_RED_COL = make_Colour_0_1( 1.0, 0.509804, 0.670588 );
        PAPAYA_WHIP_COL = make_Colour_0_1( 1.0, 0.937255, 0.835294 );
        PEACH_PUFF_COL = make_Colour_0_1( 1.0, 0.854902, 0.72549 );
        PERU_COL = make_Colour_0_1( 0.803922, 0.521569, 0.247059 );
        PINK_COL = make_Colour_0_1( 1.0, 0.752941, 0.796078 );
        PLUM_COL = make_Colour_0_1( 1.0, 0.733333, 1.0 );
        POWDER_BLUE_COL = make_Colour_0_1( 0.690196, 0.878431, 0.901961 );
        PURPLE_COL = make_Colour_0_1( 0.627451, 0.12549, 0.941176 );
        RED_COL = make_Colour_0_1( 1.0, 0.0, 0.0 );
        ROSY_BROWN_COL = make_Colour_0_1( 1.0, 0.756863, 0.756863 );
        ROYAL_BLUE_COL = make_Colour_0_1( 0.282353, 0.462745, 1.0 );
        SADDLE_BROWN_COL = make_Colour_0_1( 0.545098, 0.270588, 0.0745098 );
        SALMON_COL = make_Colour_0_1( 1.0, 0.54902, 0.411765 );
        SANDY_BROWN_COL = make_Colour_0_1( 0.956863, 0.643137, 0.376471 );
        SEASHELL_COL = make_Colour_0_1( 1.0, 0.960784, 0.933333 );
        SEA_GREEN_COL = make_Colour_0_1( 0.329412, 1.0, 0.623529 );
        SIENNA_COL = make_Colour_0_1( 1.0, 0.509804, 0.278431 );
        SKY_BLUE_COL = make_Colour_0_1( 0.529412, 0.807843, 1.0 );
        SLATE_BLUE_COL = make_Colour_0_1( 0.513725, 0.435294, 1.0 );
        SLATE_GRAY_COL = make_Colour_0_1( 0.776471, 0.886275, 1.0 );
        SLATE_GREY_COL = make_Colour_0_1( 0.439216, 0.501961, 0.564706 );
        SNOW_COL = make_Colour_0_1( 1.0, 0.980392, 0.980392 );
        SPRING_GREEN_COL = make_Colour_0_1( 0.0, 1.0, 0.498039 );
        STEEL_BLUE_COL = make_Colour_0_1( 0.388235, 0.721569, 1.0 );
        TAN_COL = make_Colour_0_1( 1.0, 0.647059, 0.309804 );
        THISTLE_COL = make_Colour_0_1( 1.0, 0.882353, 1.0 );
        TOMATO_COL = make_Colour_0_1( 1.0, 0.388235, 0.278431 );
        TURQUOISE_COL = make_Colour_0_1( 0.25098, 0.878431, 0.815686 );
        VIOLET_COL = make_Colour_0_1( 0.933333, 0.509804, 0.933333 );
        VIOLET_RED_COL = make_Colour_0_1( 1.0, 0.243137, 0.588235 );
        WHEAT_COL = make_Colour_0_1( 1.0, 0.905882, 0.729412 );
        WHITE_COL = make_Colour_0_1( 1.0, 1.0, 1.0 );
        WHITE_SMOKE_COL = make_Colour_0_1( 0.960784, 0.960784, 0.960784 );
        YELLOW_COL = make_Colour_0_1( 1.0, 1.0, 0.0 );
        YELLOW_GREEN_COL = make_Colour_0_1( 0.603922, 0.803922, 0.196078 );
        TRANSPARENT_COL = make_rgba_Colour( 0, 0, 0, 0 );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_default_surfprop
@INPUT      : 
@OUTPUT     : surfprop
@RETURNS    : 
@DESCRIPTION: Passes back a default surface property.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  get_default_surfprop(
    Surfprop  *surfprop )
{
    Surfprop_a(*surfprop) = 0.3f;
    Surfprop_d(*surfprop) = 0.3f;
    Surfprop_s(*surfprop) = 0.4f;
    Surfprop_se(*surfprop) = 10.0f;
    Surfprop_t(*surfprop) = 1.0f;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : lookup_colour
@INPUT      : colour_name
@OUTPUT     : col
@RETURNS    : TRUE if colour found
@DESCRIPTION: Given a colour name, such as "red" or "Green", finds the
              corresponding Colour.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  lookup_colour(
    STRING    colour_name,
    Colour    *col )
{
    BOOLEAN  found;
    int      i;

    check_initialize_colours();

    found = FALSE;

    for_less( i, 0, SIZEOF_STATIC_ARRAY(colour_lookup) )
    {
        if( strings_equivalent( colour_lookup[i].name, colour_name ) )
        {
            *col = *colour_lookup[i].colour;
            found = TRUE;
            break;
        }
    }

    return( found );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : lookup_colour_name
@INPUT      : col
@OUTPUT     : colour_name
@RETURNS    : TRUE if found
@DESCRIPTION: Given a colour, tries to find a named colour with the same
              values, and passes back the colour name.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  lookup_colour_name(
    Colour  col,
    STRING  *colour_name )
{
    BOOLEAN  found;
    int      i;

    check_initialize_colours();

    found = FALSE;

    for_less( i, 0, SIZEOF_STATIC_ARRAY(colour_lookup) )
    {
        if( col == *colour_lookup[i].colour )
        {
            *colour_name = create_string( colour_lookup[i].name );
            found = TRUE;
            break;
        }
    }

    return( found );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : strings_equivalent
@INPUT      : str1
              str2
@OUTPUT     : 
@RETURNS    : TRUE if strings equiv
@DESCRIPTION: Checks if two strings are equivalent, ignoring CASE and white
              space.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  BOOLEAN  strings_equivalent(
    STRING   str1,
    STRING   str2 )
{
    BOOLEAN  equivalent;
    int      i1, i2;

    equivalent = TRUE;
    i1 = 0;
    i2 = 0;

    while( str1[i1] != END_OF_STRING || str2[i2] != END_OF_STRING )
    {
        if( str1[i1] == ' ' || str1[i1] == '\t' || str1[i1] == '_' )
            ++i1;
        else if( str2[i2] == ' ' || str2[i2] == '\t' || str2[i2] == '_' )
            ++i2;
        else if( get_lower_case(str1[i1]) == get_lower_case(str2[i2]) )
        {
            ++i1;
            ++i2;
        }
        else
        {
            equivalent = FALSE;
            break;
        }
    }

    return( equivalent );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : rgb_to_hsl
@INPUT      : r
              g
              b
@OUTPUT     : h
              s
              l
@RETURNS    : 
@DESCRIPTION: Converts an RGB colour to a HSL (Hue-Saturation-Light) colour.
              This code was taken from Graphics Gems.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  rgb_to_hsl(
    Real    r,
    Real    g,
    Real    b,
    Real    *h,
    Real    *s,
    Real    *l )
{
    Real  v, m, vm, r2, g2, b2;

    check_initialize_colours();

    v = MAX3( r, g, b );

    m = MIN3( r, g, b );

    *l = (m + v) / 2.0;

    if( *l > 0.0 )
    {
        vm = v - m;
        *s = vm;

        if( *s > 0.0 )
        {
            if( *l <= 0.5 )
                *s /= v + m;
            else
                *s /= 2.0 - v - m;

            r2 = (v - r) / vm;
            g2 = (v - g) / vm;
            b2 = (v - b) / vm;

            if( r == v )
            {
                if( g == m )
                    *h = 5.0 + b2;
                else
                    *h = 1.0 - g2;
            }
            else if( g == v )
            {
                if( b == m )
                    *h = 1.0 + r2;
                else
                    *h = 3.0 - b2;
            }
            else
            {
                if( r == m )
                    *h = 3.0 + g2;
                else
                    *h = 5.0 - r2;
            }

            *h /= 6.0;
        }
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : hsl_to_rgb
@INPUT      : h
              s
              l
@OUTPUT     : r
              g
              b
@RETURNS    : 
@DESCRIPTION: Converts a HSL (Hue-Saturation-Light) colour to RGB colour.
              This code was taken from Graphics Gems.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  hsl_to_rgb(
    Real   h,
    Real   sl,
    Real   l,
    Real   *r,
    Real   *g,
    Real   *b )
{
    Real  v;
    Real  m, sv;
    int   sextant;
    Real  fract, vsf, mid1, mid2;

    if( l <= 0.5 )
    {
        v = l * (1.0 + sl);
    }
    else
    {
        v = l + sl - l * sl;
    }

    if( v <= 0.0 )
    {
        *r = 0.0;
        *g = 0.0;
        *b = 0.0;
    }
    else
    {
        m = l + l - v;
        sv = (v - m) / v;
        h *= 6.0;
        sextant = (int) h;
        fract = h - (Real) sextant;
        vsf = v * sv * fract;
        mid1 = m + vsf;
        mid2 = v - vsf;

        switch( sextant )
        {
        case 0:
        case 6:  *r = v;     *g = mid1;  *b = m;     break;
        case 1:  *r = mid2;  *g = v;     *b = m;     break;
        case 2:  *r = m;     *g = v;     *b = mid1;  break;
        case 3:  *r = m;     *g = mid2;  *b = v;     break;
        case 4:  *r = mid1;  *g = m;     *b = v;     break;
        case 5:  *r = v;     *g = m;     *b = mid2;  break;
        }
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : convert_colour_to_hsl
@INPUT      : rgb
@OUTPUT     : hsl
@RETURNS    : 
@DESCRIPTION: Converts an rgb colour to hsl.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  convert_colour_to_hsl(
    Colour   rgb,
    Colour   *hsl )
{
    Real   h, s, l;

    rgb_to_hsl( get_Colour_r_0_1(rgb),
                get_Colour_g_0_1(rgb),
                get_Colour_b_0_1(rgb),
                &h, &s, &l );

    *hsl = make_Colour_0_1( h, s, l );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : convert_colour_to_rgb
@INPUT      : hsl
@OUTPUT     : rgb
@RETURNS    : 
@DESCRIPTION: Converts a hsl colour to rgb.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  convert_colour_to_rgb(
    Colour   hsl,
    Colour   *rgb )
{
    Real   r, g, b;

    hsl_to_rgb( get_Colour_r_0_1(hsl),
                get_Colour_g_0_1(hsl),
                get_Colour_b_0_1(hsl),
                &r, &g, &b );

    *rgb = make_Colour_0_1( r, g, b );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : equal_colours
@INPUT      : col1
              col2
@OUTPUT     : 
@RETURNS    : TRUE if equal
@DESCRIPTION: Determines if two colours are equal.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  equal_colours(
    Colour  col1,
    Colour  col2 )
{
    return( col1 == col2 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : convert_colour_to_string
@INPUT      : col
@OUTPUT     : string
@RETURNS    : 
@DESCRIPTION: Converts a colour to a string.  If a named colour exists for
              the colour, then the name is assigned to the string.  Otherwise
              the 3 numerical components of the colour are assigned to the
              string, e.g. "0.3 0.7 1.0".
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  STRING  convert_colour_to_string(
    Colour   col )
{
    char       buffer[EXTREMELY_LARGE_STRING_SIZE];
    STRING     string;

    if( !lookup_colour_name( col, &string ) )
    {
        (void) sprintf( buffer, "%g %g %g",
                        get_Colour_r_0_1(col),
                        get_Colour_g_0_1(col),
                        get_Colour_b_0_1(col) );
        string = create_string( buffer );
    }

    return( string );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : convert_string_to_colour
@INPUT      : string
@OUTPUT     : 
@RETURNS    : colour
@DESCRIPTION: Converts a string to a colour.  If the string is a named colour
              such as "RED", then the corresponding colour is returned.
              Otherwise, the string is assumed to be of the form "0.3 1.0 0.7"
              and the corresponding colour returned.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Colour  convert_string_to_colour(
    STRING     string )
{
    Colour   colour;
    double   r, g, b, a;

    if( !lookup_colour( string, &colour ) )
    {
        if( sscanf( string, "%lf %lf %lf %lf", &r, &g, &b, &a ) == 4 )
            colour = make_rgba_Colour_0_1( r, g, b, a );
        else if( sscanf( string, "%lf,%lf,%lf,%lf", &r, &g, &b, &a ) == 4 )
            colour = make_rgba_Colour_0_1( r, g, b, a );
        else if( sscanf( string, "%lf %lf %lf", &r, &g, &b ) == 3 )
            colour = make_Colour_0_1( r, g, b );
        else if( sscanf( string, "%lf,%lf,%lf", &r, &g, &b ) == 3 )
            colour = make_Colour_0_1( r, g, b );
        else
            colour = make_Colour( 0, 0, 0 );
    }

    return( colour );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_colour_distance
@INPUT      : r   : in range 0 to 255
              g
              b
@OUTPUT     : 
@RETURNS    : Squared Distance
@DESCRIPTION: Returns the squared distance between two colours, one
              specified by its rgb components, and the other as a Colour
              entity.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  int  get_colour_distance(
    int      r,
    int      g,
    int      b,
    Colour   c2 )
{
    int  dr, dg, db;

    dr = r - get_Colour_r( c2 );
    dg = g - get_Colour_g( c2 );
    db = b - get_Colour_b( c2 );

    return( dr * dr + dg * dg + db * db );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : find_closest_colour
@INPUT      : r
              g
              b
              n_colours
              colours[]
@OUTPUT     : 
@RETURNS    : index into colours[]
@DESCRIPTION: Finds the closest colour in a list of colours, returning the
              index.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  int   find_closest_colour(
    int     r,
    int     g,
    int     b,
    int     n_colours,
    Colour  colours[] )
{
    int   i, dist, min_dist, min_index;

    min_index = -1;
    min_dist = 0;

    for_less( i, 0, n_colours )
    {
        dist = get_colour_distance( r, g, b, colours[i] );
        if( i == 0 || dist < min_dist )
        {
            min_index = i;
            min_dist = dist;
        }
    }

    return( min_index );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_Colour_luminance
@INPUT      : colour
@OUTPUT     : 
@RETURNS    : value between 0 and 255
@DESCRIPTION: Converts a colour to gray-scale.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  int  get_Colour_luminance(
    Colour   colour )
{
    return( ROUND( 0.299 * (Real) get_Colour_r(colour) +
                   0.587 * (Real) get_Colour_g(colour) +
                   0.114 * (Real) get_Colour_b(colour) ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : scale_colour
@INPUT      : colour
              factor
@OUTPUT     : 
@RETURNS    : Colour
@DESCRIPTION: Scales a colour by a factor, but not the alpha component.  Does
              not check for overflow.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Colour  scale_colour(
    Colour   colour,
    Real     factor )
{
    return( make_rgba_Colour( ROUND( (Real) get_Colour_r(colour) * factor ),
                              ROUND( (Real) get_Colour_g(colour) * factor ),
                              ROUND( (Real) get_Colour_b(colour) * factor ),
                              get_Colour_a(colour) ) );
}
