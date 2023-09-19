/* +++Date last modified: 05-Jul-1997 */
/* Updated comments, 05-Aug-2013 */

/*

SUNRISET.C - computes Sun rise/set times, start/end of twilight, and
             the length of the day at any date and latitude

Written as DAYLEN.C, 1989-08-16

Modified to SUNRISET.C, 1992-12-01

(c) Paul Schlyter, 1989, 1992

Released to the public domain by Paul Schlyter, December 1992

*/


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <getopt.h>
#include <math.h>


/* A macro to compute the number of days elapsed since 2000 Jan 0.0 */
/* (which is equal to 1999 Dec 31, 0h UT)                           */

#define days_since_2000_Jan_0(y,m,d) \
    (367L*(y)-((7*((y)+(((m)+9)/12)))/4)+((275*(m))/9)+(d)-730530L)

/* Some conversion factors between radians and degrees */

#ifndef PI
 #define PI        3.1415926535897932384
#endif

#define RADEG     ( 180.0 / PI )
#define DEGRAD    ( PI / 180.0 )

/* The trigonometric functions in degrees */

#define sind(x)  sin((x)*DEGRAD)
#define cosd(x)  cos((x)*DEGRAD)
#define tand(x)  tan((x)*DEGRAD)

#define atand(x)    (RADEG*atan(x))
#define asind(x)    (RADEG*asin(x))
#define acosd(x)    (RADEG*acos(x))
#define atan2d(y,x) (RADEG*atan2(y,x))


/* Following are some macros around the "workhorse" function __daylen__ */
/* They mainly fill in the desired values for the reference altitude    */
/* below the horizon, and also selects whether this altitude should     */
/* refer to the Sun's center or its upper limb.                         */


/* This macro computes the length of the day, from sunrise to sunset. */
/* Sunrise/set is considered to occur when the Sun's upper limb is    */
/* 35 arc minutes below the horizon (this accounts for the refraction */
/* of the Earth's atmosphere).                                        */
#define day_length(year,month,day,lon,lat)  \
        __daylen__( year, month, day, lon, lat, -35.0/60.0, 1 )

/* This macro computes the length of the day, including civil twilight. */
/* Civil twilight starts/ends when the Sun's center is 6 degrees below  */
/* the horizon.                                                         */
#define day_civil_twilight_length(year,month,day,lon,lat)  \
        __daylen__( year, month, day, lon, lat, -6.0, 0 )

/* This macro computes the length of the day, incl. nautical twilight.  */
/* Nautical twilight starts/ends when the Sun's center is 12 degrees    */
/* below the horizon.                                                   */
#define day_nautical_twilight_length(year,month,day,lon,lat)  \
        __daylen__( year, month, day, lon, lat, -12.0, 0 )

/* This macro computes the length of the day, incl. astronomical twilight. */
/* Astronomical twilight starts/ends when the Sun's center is 18 degrees   */
/* below the horizon.                                                      */
#define day_astronomical_twilight_length(year,month,day,lon,lat)  \
        __daylen__( year, month, day, lon, lat, -18.0, 0 )


/* This macro computes times for sunrise/sunset.                      */
/* Sunrise/set is considered to occur when the Sun's upper limb is    */
/* 35 arc minutes below the horizon (this accounts for the refraction */
/* of the Earth's atmosphere).                                        */
#define sun_rise_set(year,month,day,lon,lat,tz,rise,set)  \
        __sunriset__( year, month, day, lon, lat, tz, -35.0/60.0, 1, rise, set )

/* This macro computes the start and end times of civil twilight.       */
/* Civil twilight starts/ends when the Sun's center is 6 degrees below  */
/* the horizon.                                                         */
#define civil_twilight(year,month,day,lon,lat,tz,start,end)  \
        __sunriset__( year, month, day, lon, lat, tz, -6.0, 0, start, end )

/* This macro computes the start and end times of nautical twilight.    */
/* Nautical twilight starts/ends when the Sun's center is 12 degrees    */
/* below the horizon.                                                   */
#define nautical_twilight(year,month,day,lon,lat,tz,start,end)  \
        __sunriset__( year, month, day, lon, lat, tz, -12.0, 0, start, end )

/* This macro computes the start and end times of astronomical twilight.   */
/* Astronomical twilight starts/ends when the Sun's center is 18 degrees   */
/* below the horizon.                                                      */
#define astronomical_twilight(year,month,day,lon,lat,tz,start,end)  \
        __sunriset__( year, month, day, lon, lat, tz, -18.0, 0, start, end )


/* Function prototypes */

double __daylen__( int year, int month, int day, double lon, double lat,
                   double altit, int upper_limb );

int __sunriset__( int year, int month, int day, double lon, double lat, double tz,
                  double altit, int upper_limb, double *rise, double *set );

void sunpos( double d, double *lon, double *r );

void sun_RA_dec( double d, double *RA, double *dec, double *r );

double revolution( double x );

double rev180( double x );

double GMST0( double d );

int H( double r )
{
  return (int)r;
}

int M( double r )
{
  return (int)((r - (int)r)*60);
}



/* A small test program */

static void print_usage(const char *prog)
{
        printf("Usage: %s [-ymdtoah]\n", prog);
        puts("  -y --year     year (default today)\n"
             "  -m --mon      month (default today)\n"
             "  -d --day      day (default today)\n"
             "  -t --tz       time zone(default 9)\n"
             "  -o --lon      longitude (default 126.978420)\n"
             "  -a --lat      latitude (default 37.566692) \n"
             "  -h --help     this message\n");
        exit(1);
}

static void parse_opts(int argc, char *argv[], int *y, int *m, int *d, double *lon, double *lat, double *tz)
{
  time_t systime;
  time(&systime);
  struct tm * now = localtime(&systime);

  *y = now->tm_year + 1900;
  *m = now->tm_mon + 1;
  *d = now->tm_mday;
  *lon = 126.978420;
  *lat = 37.566692;
  *tz = 9;

  while(1)
  {
    static const struct option lopts[] = {
      { "year", required_argument, 0, 'y' },
      { "mon",  required_argument, 0, 'm' },
      { "day",  required_argument, 0, 'd' },
      { "tz",   required_argument, 0, 't' },
      { "lon",  required_argument, 0, 'o' },
      { "lat",  required_argument, 0, 'a' },
      { "help", no_argument,       0, 'h' },
      { NULL,   0,                 0,  0 }
    };

    int c;

    c = getopt_long(argc, argv, "y:m:d:t:o:a:", lopts, NULL);
    if(c == -1) break;

    switch(c)
    {
      case 'y': sscanf(optarg, "%d", y); break;
      case 'm': sscanf(optarg, "%d", m); break;
      case 'd': sscanf(optarg, "%d", d); break;
      case 'o': sscanf(optarg, "%lf", lon); break;
      case 'a': sscanf(optarg, "%lf", lat); break;
      case 't': sscanf(optarg, "%lf", tz); break;
      case 'h':
      default: print_usage(argv[0]); break;
    }
  }
}

int main(int argc, char * argv[])
{
	int year,month,day;
	double lon, lat;
	double tz;
	double daylen, civlen, nautlen, astrlen;
	double rise, set, civ_start, civ_end, naut_start, naut_end, astr_start, astr_end;
	int rs, civ, naut, astr;

	parse_opts(argc, argv, &year, &month, &day, &lon, &lat, &tz);

	printf("Longitude:       %c%.6lf\n", (lon<0)?'W':'E', (lon<0)?-lon:lon);
	printf("Latitude:        %c%.6lf\n", (lat<0)?'S':'N', (lat<0)?-lat:lat);
	printf("Date:            %d %02d %02d\n", year, month, day);
	printf("Timezone:        %d\n", (int)tz);

	daylen  = day_length(year,month,day,lon,lat);
	civlen  = day_civil_twilight_length(year,month,day,lon,lat);
	nautlen = day_nautical_twilight_length(year,month,day,lon,lat);
	astrlen = day_astronomical_twilight_length(year,month,day,lon,lat);

	printf( "With civil twilight         %02d:%02d\n", H(civlen), M(civlen));
	printf( "Day length:                 %02d:%02d\n", H(daylen), M(daylen));
	printf( "With nautical twilight      %02d:%02d\n", H(nautlen), M(nautlen));
	printf( "With astronomical twilight  %02d:%02d\n", H(astrlen), M(astrlen));
	printf( "Length of twilight: civil   %02d:%02d\n", H((civlen-daylen)/2.0), M((civlen-daylen)/2.0));
	printf( "                 nautical   %02d:%02d\n", H((nautlen-daylen)/2.0), M((nautlen-daylen)/2.0));
	printf( "             astronomical   %02d:%02d\n", H((astrlen-daylen)/2.0), M((astrlen-daylen)/2.0));

	rs   = sun_rise_set         ( year, month, day, lon, lat, tz, &rise, &set );
	civ  = civil_twilight       ( year, month, day, lon, lat, tz, &civ_start, &civ_end );
	naut = nautical_twilight    ( year, month, day, lon, lat, tz, &naut_start, &naut_end );
	astr = astronomical_twilight( year, month, day, lon, lat, tz, &astr_start, &astr_end );

	printf( "Sun at south %02d:%02d\n", H((rise+set)/2.0), M((rise+set)/2.0));

	switch( rs )
	{
		case 0:
			printf( "Sun rises %02d:%02d, sets %02d:%02d\n", H(rise+tz), M(rise+tz), H(set+tz), M(set+tz));
			break;
		case +1:
			printf( "Sun above horizon\n" );
			break;
		case -1:
			printf( "Sun below horizon\n" );
			break;
	}

	switch( civ )
	{
		case 0:
			printf( "Civil twilight starts %02d:%02d, ends %02d:%02d\n", H(civ_start+tz), M(civ_start+tz), H(civ_end+tz), M(civ_end+tz));
			break;
		case +1:
			printf( "Never darker than civil twilight\n" );
			break;
		case -1:
			printf( "Never as bright as civil twilight\n" );
			break;
	}

	switch( naut )
	{
		case 0:
			printf( "Nautical twilight starts %02d:%02d, ends %02d:%02d\n", H(naut_start+tz), M(naut_start+tz), H(naut_end+tz), M(naut_end+tz));
			break;
		case +1:
			printf( "Never darker than nautical twilight\n" );
			break;
		case -1:
			printf( "Never as bright as nautical twilight\n" );
			break;
	}

	switch( astr )
	{
		case 0:
			printf( "Astronomical twilight starts %02d:%02d, ends %02d:%02d\n", H(astr_start+tz), M(astr_start+tz), H(astr_end+tz), M(astr_end+tz));
			break;
		case +1:
			printf( "Never darker than astronomical twilight\n" );
			break;
		case -1:
			printf( "Never as bright as astronomical twilight\n" );
			break;
	}
}


/* The "workhorse" function for sun rise/set times */

int __sunriset__( int year, int month, int day, double lon, double lat, double tz,
                  double altit, int upper_limb, double *trise, double *tset )
/***************************************************************************/
/* Note: year,month,date = calendar date, 1801-2099 only.             */
/*       Eastern longitude positive, Western longitude negative       */
/*       Northern latitude positive, Southern latitude negative       */
/*       The longitude value IS critical in this function!            */
/*       altit = the altitude which the Sun should cross              */
/*               Set to -35/60 degrees for rise/set, -6 degrees       */
/*               for civil, -12 degrees for nautical and -18          */
/*               degrees for astronomical twilight.                   */
/*         upper_limb: non-zero -> upper limb, zero -> center         */
/*               Set to non-zero (e.g. 1) when computing rise/set     */
/*               times, and to zero when computing start/end of       */
/*               twilight.                                            */
/*        *rise = where to store the rise time                        */
/*        *set  = where to store the set  time                        */
/*                Both times are relative to the specified altitude,  */
/*                and thus this function can be used to compute       */
/*                various twilight times, as well as rise/set times   */
/* Return value:  0 = sun rises/sets this day, times stored at        */
/*                    *trise and *tset.                               */
/*               +1 = sun above the specified "horizon" 24 hours.     */
/*                    *trise set to time when the sun is at south,    */
/*                    minus 12 hours while *tset is set to the south  */
/*                    time plus 12 hours. "Day" length = 24 hours     */
/*               -1 = sun is below the specified "horizon" 24 hours   */
/*                    "Day" length = 0 hours, *trise and *tset are    */
/*                    both set to the time when the sun is at south.  */
/*                                                                    */
/**********************************************************************/
{
	double  d,  /* Days since 2000 Jan 0.0 (negative before) */
	sr,         /* Solar distance, astronomical units */
	sRA,        /* Sun's Right Ascension */
	sdec,       /* Sun's declination */
	sradius,    /* Sun's apparent radius */
	t,          /* Diurnal arc */
	tsouth,     /* Time when Sun is at south */
	sidtime;    /* Local sidereal time */

	int rc = 0; /* Return cde from function - usually 0 */

	/* Compute d of 12h local mean solar time */
	d = days_since_2000_Jan_0(year,month,day) + 0.5 - lon/360.0;

	/* Compute the local sidereal time of this moment */
	sidtime = revolution( GMST0(d) + 180.0 + lon );

	/* Compute Sun's RA, Decl and distance at this moment */
	sun_RA_dec( d, &sRA, &sdec, &sr );

	/* Compute time when Sun is at south - in hours UT */
	tsouth = 12.0 - rev180(sidtime - sRA)/15.0;

	/* Compute the Sun's apparent radius in degrees */
	sradius = 0.2666 / sr;

	/* Do correction to upper limb, if necessary */
	if ( upper_limb ) altit -= sradius;

	/* Compute the diurnal arc that the Sun traverses to reach */
	/* the specified altitude altit: */
	{
		double cost;
		cost = ( sind(altit) - sind(lat) * sind(sdec) ) / ( cosd(lat) * cosd(sdec) );
		if ( cost >= 1.0 )
			rc = -1, t = 0.0;       /* Sun always below altit */
		else if ( cost <= -1.0 )
			rc = +1, t = 12.0;      /* Sun always above altit */
		else
			t = acosd(cost)/15.0;   /* The diurnal arc, hours */
	}

	/* Store rise and set times - in hours UT */
	*trise = tsouth - t;
	*tset  = tsouth + t;

	return rc;
}  /* __sunriset__ */



/* The "workhorse" function */


double __daylen__( int year, int month, int day, double lon, double lat,
                   double altit, int upper_limb )
/**********************************************************************/
/* Note: year,month,date = calendar date, 1801-2099 only.             */
/*       Eastern longitude positive, Western longitude negative       */
/*       Northern latitude positive, Southern latitude negative       */
/*       The longitude value is not critical. Set it to the correct   */
/*       longitude if you're picky, otherwise set to to, say, 0.0     */
/*       The latitude however IS critical - be sure to get it correct */
/*       altit = the altitude which the Sun should cross              */
/*               Set to -35/60 degrees for rise/set, -6 degrees       */
/*               for civil, -12 degrees for nautical and -18          */
/*               degrees for astronomical twilight.                   */
/*         upper_limb: non-zero -> upper limb, zero -> center         */
/*               Set to non-zero (e.g. 1) when computing day length   */
/*               and to zero when computing day+twilight length.      */
/**********************************************************************/
{
	double  d,  /* Days since 2000 Jan 0.0 (negative before) */
	obl_ecl,    /* Obliquity (inclination) of Earth's axis */
	sr,         /* Solar distance, astronomical units */
	slon,       /* True solar longitude */
	sin_sdecl,  /* Sine of Sun's declination */
	cos_sdecl,  /* Cosine of Sun's declination */
	sradius;    /* Sun's apparent radius */
	double t;     /* Diurnal arc */

	/* Compute d of 12h local mean solar time */
	d = days_since_2000_Jan_0(year,month,day) + 0.5 - lon/360.0;

	/* Compute obliquity of ecliptic (inclination of Earth's axis) */
	obl_ecl = 23.4393 - 3.563E-7 * d;

	/* Compute Sun's ecliptic longitude and distance */
	sunpos( d, &slon, &sr );

	/* Compute sine and cosine of Sun's declination */
	sin_sdecl = sind(obl_ecl) * sind(slon);
	cos_sdecl = sqrt( 1.0 - sin_sdecl * sin_sdecl );

	/* Compute the Sun's apparent radius, degrees */
	sradius = 0.2666 / sr;

	/* Do correction to upper limb, if necessary */
	if ( upper_limb ) altit -= sradius;

	/* Compute the diurnal arc that the Sun traverses to reach */
	/* the specified altitude altit: */
	{
		double cost;
		cost = ( sind(altit) - sind(lat) * sin_sdecl ) /
			( cosd(lat) * cos_sdecl );
		if ( cost >= 1.0 )
			t = 0.0;                      /* Sun always below altit */
		else if ( cost <= -1.0 )
			t = 24.0;                     /* Sun always above altit */
		else  t = (2.0/15.0) * acosd(cost); /* The diurnal arc, hours */
	}

	return t;
}  /* __daylen__ */


/* This function computes the Sun's position at any instant */

void sunpos( double d, double *lon, double *r )
/******************************************************/
/* Computes the Sun's ecliptic longitude and distance */
/* at an instant given in d, number of days since     */
/* 2000 Jan 0.0.  The Sun's ecliptic latitude is not  */
/* computed, since it's always very near 0.           */
/******************************************************/
{
	double M,         /* Mean anomaly of the Sun */
		w,         /* Mean longitude of perihelion */
					/* Note: Sun's mean longitude = M + w */
		e,         /* Eccentricity of Earth's orbit */
		E,         /* Eccentric anomaly */
		x, y,      /* x, y coordinates in orbit */
		v;         /* True anomaly */

	/* Compute mean elements */
	M = revolution( 356.0470 + 0.9856002585 * d );
	w = 282.9404 + 4.70935E-5 * d;
	e = 0.016709 - 1.151E-9 * d;

	/* Compute true longitude and radius vector */
	E = M + e * RADEG * sind(M) * ( 1.0 + e * cosd(M) );
	x = cosd(E) - e;
	y = sqrt( 1.0 - e*e ) * sind(E);
	*r = sqrt( x*x + y*y );              /* Solar distance */
	v = atan2d( y, x );                  /* True anomaly */
	*lon = v + w;                        /* True solar longitude */
	if ( *lon >= 360.0 )
		*lon -= 360.0;                   /* Make it 0..360 degrees */
}

void sun_RA_dec( double d, double *RA, double *dec, double *r )
/******************************************************/
/* Computes the Sun's equatorial coordinates RA, Decl */
/* and also its distance, at an instant given in d,   */
/* the number of days since 2000 Jan 0.0.             */
/******************************************************/
{
	double lon, obl_ecl, x, y, z;

	/* Compute Sun's ecliptical coordinates */
	sunpos( d, &lon, r );

	/* Compute ecliptic rectangular coordinates (z=0) */
	x = *r * cosd(lon);
	y = *r * sind(lon);

	/* Compute obliquity of ecliptic (inclination of Earth's axis) */
	obl_ecl = 23.4393 - 3.563E-7 * d;

	/* Convert to equatorial rectangular coordinates - x is unchanged */
	z = y * sind(obl_ecl);
	y = y * cosd(obl_ecl);

	/* Convert to spherical coordinates */
	*RA = atan2d( y, x );
	*dec = atan2d( z, sqrt(x*x + y*y) );

}  /* sun_RA_dec */


/******************************************************************/
/* This function reduces any angle to within the first revolution */
/* by subtracting or adding even multiples of 360.0 until the     */
/* result is >= 0.0 and < 360.0                                   */
/******************************************************************/

#define INV360    ( 1.0 / 360.0 )

double revolution( double x )
/*****************************************/
/* Reduce angle to within 0..360 degrees */
/*****************************************/
{
	return( x - 360.0 * floor( x * INV360 ) );
}  /* revolution */

double rev180( double x )
/*********************************************/
/* Reduce angle to within +180..+180 degrees */
/*********************************************/
{
	return( x - 360.0 * floor( x * INV360 + 0.5 ) );
}  /* revolution */


/*******************************************************************/
/* This function computes GMST0, the Greenwich Mean Sidereal Time  */
/* at 0h UT (i.e. the sidereal time at the Greenwhich meridian at  */
/* 0h UT).  GMST is then the sidereal time at Greenwich at any     */
/* time of the day.  I've generalized GMST0 as well, and define it */
/* as:  GMST0 = GMST - UT  --  this allows GMST0 to be computed at */
/* other times than 0h UT as well.  While this sounds somewhat     */
/* contradictory, it is very practical:  instead of computing      */
/* GMST like:                                                      */
/*                                                                 */
/*  GMST = (GMST0) + UT * (366.2422/365.2422)                      */
/*                                                                 */
/* where (GMST0) is the GMST last time UT was 0 hours, one simply  */
/* computes:                                                       */
/*                                                                 */
/*  GMST = GMST0 + UT                                              */
/*                                                                 */
/* where GMST0 is the GMST "at 0h UT" but at the current moment!   */
/* Defined in this way, GMST0 will increase with about 4 min a     */
/* day.  It also happens that GMST0 (in degrees, 1 hr = 15 degr)   */
/* is equal to the Sun's mean longitude plus/minus 180 degrees!    */
/* (if we neglect aberration, which amounts to 20 seconds of arc   */
/* or 1.33 seconds of time)                                        */
/*                                                                 */
/*******************************************************************/

double GMST0( double d )
{
	double sidtim0;
	/* Sidtime at 0h UT = L (Sun's mean longitude) + 180.0 degr  */
	/* L = M + w, as defined in sunpos().  Since I'm too lazy to */
	/* add these numbers, I'll let the C compiler do it for me.  */
	/* Any decent C compiler will add the constants at compile   */
	/* time, imposing no runtime or code overhead.               */
	sidtim0 = revolution( ( 180.0 + 356.0470 + 282.9404 ) +
							( 0.9856002585 + 4.70935E-5 ) * d );
	return sidtim0;
}  /* GMST0 */
