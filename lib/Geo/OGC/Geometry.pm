## @namespace Geo::OGC
# @brief The simple feature geometries
#
# The classes and methods in the Geo::OGC:: namespace should conform
# to the OGC implementation specification (currently 06-103r3) for
# Geographic Information - Simple feature access.
# @note Some of the brief descriptions in this document are taken
# verbatim from <a
# href="http://www.opengeospatial.org/standards/sfa">the OGC
# specification (OGC 06-103r3 Version 1.2.0)</a>. <a
# href="http://www.opengeospatial.org/ogc/document">Copyright (c) 2006
# Open Geospatial Consortium, Inc. All Rights Reserved.</a>
# @note Calling a method that is not yet implemented causes an
# exception. The status of implementation is not always shown in this
# document.
# @note Most of the methods for testing spatial relations and the
# methods that support spatial analysis are not yet implemented.

package Geo::OGC::Geometry;

=pod

=head1 NAME

Geo::OGC::Geometry - Simple feature geometry classes

The <a
href="http://map.hut.fi/doc/Geoinformatica/html/class_geo_1_1_o_g_c_1_1_geometry.html">
documentation of Geo::OGC::Geometry</a> is written in doxygen format.

=cut

use strict;
use Carp;

BEGIN {
    use Exporter "import";
    our @ISA = qw(Exporter);
    our %EXPORT_TAGS = ( 'all' => [ qw( &ccw &intersect &distance_point_line ) ] );
    our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
    our $VERSION = '0.01';
}

## @cmethod Geo::OGC::Geometry new(%params)
# @brief Create a new and initialized object.
#
# @param params A set of named parameters (the subclasses may add new
# known parameters):
# - <i>Text</i> A well-known text, constructs an object of class the
# text defines.
# - <i>SRID</i> Sets the SRID attribute of the object, default is -1.
# - <i>Precision</i> If specified, has the effect that numeric
# comparisons in the Equals method is is preceded with a rounding
# operation (using sprintf "%.pe", where p is the Precision-1, i.e the
# number of meaningful numbers is Precision). Affects also AsText.
#
# This method should be called eventually by all constructors. Blesses
# the object into the correct class and calls init.
sub new {
    my($package, %params) = @_;
    return parse_wkt($params{Text}) if exists $params{Text};
    my $self = {};
    bless $self => (ref($package) or $package);
    $self->init(%params);
    return $self;
}

## @method void init(%params)
# @brief Override in new classes but call $self->SUPER::init(%params);
sub init {
    my($self, %params) = @_;
    $self->{SRID} = exists $params{SRID} ? $params{SRID} : -1;
    $self->{Precision} = $params{Precision}-1 if exists $params{Precision};
}

## @method void copy($clone)
# @brief Copies the attributes of this object to the other object,
# which is a newly created object of same class.
sub copy {
    my($self, $clone) = @_;
    $clone->{SRID} = $self->{SRID};
    $clone->{Precision} = $self->{Precision} if exists $self->{Precision};
}

## @method Geo::OGC::Geometry parse_wkt($Text)
# @brief parse well known text and construct respective geometry
sub parse_wkt {
    my $self;
    for ($_[0]) {
	if (/^\s*POINT/i) {
	    s/^\s*POINT\s*([ZM]*)\s*\(\s*//i;
	    my $m = lc($1);
	    s/\s*\)\s*$//;
	    my @point = split /\s+/;
	    $self = Geo::OGC::Point->new("point$m"=>\@point);
	} elsif (/^\s*MULTIPOINT/i) {
	    s/^\s*MULTIPOINT\s*([ZM]*)\s*\(\s*//i;
	    my $m = $1;
	    s/\s*\)\s*$//;
	    my @points = split /\s*,\s*/;
	    $self = Geo::OGC::MultiPoint->new();
	    for (@points) {
		$self->AddGeometry(parse_wkt("POINT$m ($_)"));
	    }
	} elsif (/^\s*LINESTRING/i) {
	    s/^\s*LINESTRING\s*([ZM]*)\s*\(\s*//i;
	    my $m = $1;
	    s/\s*\)\s*$//;
	    my @points = split /\s*,\s*/;
	    $self = Geo::OGC::LineString->new();
	    for (@points) {
		$self->AddPoint(parse_wkt("POINT$m ($_)"));
	    }
	} elsif (/^\s*MULTILINESTRING/i) {
	    s/^\s*MULTILINESTRING\s*([ZM]*)[\s\(]*//i;
	    my $m = $1;
	    s/[\s\)]*$//;
	    my @strings = split /\)\s*,\s*\(/;
	    $self = Geo::OGC::MultiLineString->new();
	    for (@strings) {
		$self->AddGeometry(parse_wkt("LINESTRING$m ($_)"));
	    }
	} elsif (/^\s*LINEARRING/i) {
	    s/^\s*LINEARRING\s*([ZM]*)\s*\(\s*//i;
	    my $m = $1;
	    s/\s*\)\s*$//;
	    my @points = split /\s*,\s*/;
	    $self = Geo::OGC::LinearRing->new();
	    for (@points) {
		$self->AddPoint(parse_wkt("POINT$m ($_)"));
	    }
	} elsif (/^\s*POLYGON/i) {
	    s/^\s*POLYGON\s*([ZM]*)[\s\(]*//i;
	    my $m = $1;
	    s/[\s\)]*$//;
	    my @rings = split /\)\s*,\s*\(/;
	    $self = Geo::OGC::Polygon->new();
	    $self->ExteriorRing(parse_wkt("LINEARRING$m (".shift(@rings).")"));
	    for (@rings) {
		$self->AddInteriorRing(parse_wkt("LINEARRING$m ($_)"));
	    }
	} elsif (/^\s*POLYHEDRALSURFACE/i) {
	    s/^\s*POLYHEDRALSURFACE\s*([ZM]*)[\s\(]*//i;
	    my $m = $1;
	    s/[\s\)]*$//;
	    my @patches = split /\)\s*,\s*\(/;
	    $self = Geo::OGC::PolyhedralSurface->new();
	    for (@patches) {
		$self->AddPatch(parse_wkt("LINEARRING$m ($_)"));
	    }
	} elsif (/^\s*MULTIPOLYGON/i) {
	    s/^\s*MULTIPOLYGON\s*([ZM]*)[\s\(]*//i;
	    my $m = $1;
	    s/[\s\)]*$//;
	    my @polygons = split /\)\s*\)\s*,\s*\(\s*\(/;
	    $self = Geo::OGC::MultiPolygon->new();
	    for (@polygons) {
		$self->AddGeometry(parse_wkt("POLYGON$m (($_))"));
	    }
	} elsif (/^\s*GEOMETRYCOLLECTION/i) {
	    s/^\s*GEOMETRYCOLLECTION\s*([ZM]*)\s*\(\s*//i;
	    my $m = $1;
	    s/\s*\)\s*$//;
	    my @b = /,([A-Z])/g;
	    unshift @b,'';
	    my @geometries = split /,[A-Z]/;
	    $self = Geo::OGC::GeometryCollection->new();
	    for my $i (0..$#geometries) {
		$self->AddGeometry(parse_wkt($b[$i].$geometries[$i]));
	    }
	} else {
	    my $beginning = substr $_, 0, 20;
	    croak "can't parse text beginning as: $beginning";
	}
    }
    return $self;
}

## @fn $ccw($x0, $y0, $x1, $y1, $x2, $y2)
# @brief counterclockwise from Sedgewick: Algorithms in C
# @return -1, 0, or 1
sub ccw {
    my($x0, $y0, $x1, $y1, $x2, $y2) = @_;
    my $dx1 = $x1 - $x0; 
    my $dy1 = $y1 - $y0;
    my $dx2 = $x2 - $x0; 
    my $dy2 = $y2 - $y0;
    return +1 if $dx1*$dy2 > $dy1*$dx2;
    return -1 if $dx1*$dy2 < $dy1*$dx2;
    return -1 if ($dx1*$dx2 < 0) or ($dy1*$dy2 < 0);
    return +1 if ($dx1*$dx1+$dy1*$dy1) < ($dx2*$dx2+$dy2*$dy2);
    return 0;
}

## @fn $intersect($x11, $y11, $x12, $y12, $x21, $y21, $x22, $y22)
# @brief Test intersection of two lines from Sedgewick: Algorithms in C
sub intersect { # also from Sedgewick: Algorithms in C
    my($x11, $y11, $x12, $y12, $x21, $y21, $x22, $y22) = @_;
    return ((ccw($x11, $y11, $x12, $y12, $x21, $y21)
	     *ccw($x11, $y11, $x12, $y12, $x22, $y22)) <= 0)
	&& ((ccw($x21, $y21, $x22, $y22, $x11, $y11)
	     *ccw($x21, $y21, $x22, $y22, $x12, $y12)) <= 0);
}

## @fn $distance_point_line($x, $y, $x1, $y1, $x2, $y2)
# @brief Compute the distance of a point to a line.
sub distance_point_line {
    my($x, $y, $x1, $y1, $x2, $y2) = @_;
    my $l2 = ($x2-$x1)**2 + ($y2-$y1)**2;
    my $u = (($x - $x1) * ($x2 - $x1) + ($y - $y1) * ($y2 - $y1)) / $l2;
    if ($u < 0) { # distance to point 1
	return sqrt(($x-$x1)**2 + ($y-$y1)**2);
    } elsif ($u > 1) { # distance to point 2
	return sqrt(($x-$x2)**2 + ($y-$y2)**2);
    } else {
	my $ix = $x1 + $u * ($x2 - $x1);
	my $iy = $y1 + $u * ($y2 - $y1);
	return sqrt(($x-$ix)**2 + ($y-$iy)**2);
    }
}

## @ignore
sub dump {
    my($self) = @_;
    print "$self\n";
    for (sort keys %$self) {
	print "$_ => $self->{$_}\n";
    }
}

## @method Geo::OGC::Geometry Clone()
# @brief Clones this object, i.e., creates a spatially exact copy.
sub Clone {
    my($self) = @_;
    my $clone = Geo::OGC::Geometry::new($self);
    $self->copy($clone);
    return $clone;
}

## @method $Precision($Precision)
# @brief Sets or gets the precision
# @note Not in the specification
sub Precision {
    my($self, $Precision) = @_;
    defined $Precision ? 
	$self->{Precision} = $Precision-1 : $self->{Precision}+1;
}

## @method $Dimension()
# @brief The dimension of this geometric object. In non-homogeneous
# collections, this will return the largest topological dimension of
# the contained objects.
# @return 2 or 3
sub Dimension {
    my($self) = @_;
    return $self->Is3D ? 3 : 2;
}

## @method $GeometryType()
# @brief Returns the name of the class of this geometric object.
sub GeometryType {
    my($self) = @_;
    croak "GeometryType method for class ".ref($self)." is not implemented yet";
}

## @method $SRID($SRID)
# @brief Returns (or sets) the Spatial Reference System ID for this
# geometric object.
# @param SRID [optional] The new SRID if setting it.
sub SRID {
    my($self, $SRID) = @_;
    defined $SRID ? 
	$self->{SRID} = $SRID : $self->{SRID};
}

## @method Geo::OGC::Polygon Envelope()
# @brief The minimum bounding box for this Geometry.
# @note This library returns always the envelope as a ring [(minx,
# miny), (maxx, miny), (maxx, maxy), (minx, maxy), (minx, miny)]
sub Envelope {
    my($self) = @_;
    croak "Envelope method for class ".ref($self)." is not implemented yet";
}

## @method $as_text($force_parens, $include_tag)
# @brief A helper method used by AsText
sub as_text {
    my($self, $force_parens, $include_tag) = @_;
    croak "as_text method for class ".ref($self)." is not implemented yet";
}

## @method $AsText()
# @brief Returns this object in Well-known Text Representation of Geometry
sub AsText {
    my($self) = @_;
    return $self->as_text(1, 1);
}

## @method $AsBinary()
# @brief Returns this object in Well-known Binary Representation of Geometry
# @note Not implemented yet.
sub AsBinary {
    my($self) = @_;
    croak "AsBinary method for class ".ref($self)." is not implemented yet";
}

## @method $IsEmpty()
# @brief Returns true if this object is empty, i.e. contains no points
# or no data.
sub IsEmpty {
    my($self) = @_;
    croak "IsEmpty method for class ".ref($self)." is not implemented yet";
}

## @method $IsSimple()
# @brief Returns true if this geometric object has no anomalous
# geometric points, such as self intersection or self tangency.
sub IsSimple {
    my($self) = @_;
    croak "IsSimple method for class ".ref($self)." is not implemented yet";
}

## @method $Is3D()
# @brief Returns true if this geometric object has z coordinate values.
sub Is3D {
    my($self) = @_;
    croak "Is3D method for class ".ref($self)." is not implemented yet";
}

## @method $IsMeasured()
# @brief Returns true if this geometric object has m coordinate
# values.
sub IsMeasured {
    my($self) = @_;
    croak "IsMeasured method for class ".ref($self)." is not implemented yet";
}

## @method Geo::OGC::Geometry Boundary()
# @brief Returns the closure of the combinatorial boundary of this
# geometric object.
# @note Not implemented yet.
sub Boundary {
    my($self) = @_;
    croak "Boundary method for class ".ref($self)." is not implemented yet";
}

## @method $Equals($geom)
# @brief Returns true if this geometric object is "spatially equal" to
# the another geometry.
sub Equals {
    my($self, $geom) = @_;
    croak "Equals method for class ".ref($self)." is not implemented yet";
}

sub Disjoint {
    my($self, $geom) = @_;
    croak "Disjoint method for class ".ref($self)." is not implemented yet";
}

sub Intersects {
    my($self, $geom) = @_;
    croak "Intersects method for class ".ref($self)." is not implemented yet";
}

sub Touches {
    my($self, $geom) = @_;
    croak "Touches method for class ".ref($self)." is not implemented yet";
}

sub Crosses {
    my($self, $geom) = @_;
    croak "Crosses method for class ".ref($self)." is not implemented yet";
}

sub Within {
    my($self, $geom) = @_;
    croak "Within method for class ".ref($self)." is not implemented yet";
}

sub Contains {
    my($self, $geom) = @_;
    croak "Contains method for class ".ref($self)." is not implemented yet";
}

sub Overlaps {
    my($self, $geom) = @_;
    croak "Overlaps method for class ".ref($self)." is not implemented yet";
}

sub Relate {
    my($self, $geom, $int_pat) = @_;
    croak "Relate method for class ".ref($self)." is not implemented yet";
}

sub LocateAlong {
    my($self, $mValue) = @_;
    croak "LocateAlong method for class ".ref($self)." is not implemented yet";
}

sub LocateBetween {
    my($self, $mStart, $mEnd) = @_;
    croak "LocateBetween method for class ".ref($self)." is not implemented yet";
}

sub Distance {
    my($self, $geom) = @_;
    croak "Distance method for class ".ref($self)." is not implemented yet";
}

sub Buffer {
    my($self, $distance) = @_;
    croak "Buffer method for class ".ref($self)." is not implemented yet";
}

sub ConvexHull {
    my($self) = @_;
    croak "ConvexHull method for class ".ref($self)." is not implemented yet";
}

sub Intersection {
    my($self, $geom) = @_;
    croak "Intersection method for class ".ref($self)." is not implemented yet";
}

sub Union {
    my($self, $geom) = @_;
    croak "Union method for class ".ref($self)." is not implemented yet";
}

sub Difference {
    my($self, $geom) = @_;
    croak "Difference method for class ".ref($self)." is not implemented yet";
}

sub SymDifference {
    my($self, $geom) = @_;
    croak "SymDifference method for class ".ref($self)." is not implemented yet";
}

## @method MakeCollection()
# @brief Creates a collection which contains this geometry
# @note Not in the specification
sub MakeCollection {
    my($self) = @_;
    croak "MakeCollection method for class ".ref($self)." is not implemented";
}

## @method ApplyTransformation(transf)
# @param transf A point transformation method which will be applied
# for all points in the geometry as:
# @code
# ($new_x, $new_y, $new_z) = $transf->($x, $y, $z)
# @endcode
# @note Not in the specification
sub ApplyTransformation {
    my($self, $transf) = @_;
    croak "ApplyTransformation method for class ".ref($self)." is not implemented";
}

#
#    SpatialReferenceSystem
#

package Geo::OGC::SpatialReferenceSystem;

use strict;
use Carp;

sub new {
    my($package, %params) = @_;
    my $self = {};
    bless $self => (ref($package) or $package);
}

#
#    Point
#

package Geo::OGC::Point;

use strict;
use UNIVERSAL qw(isa);
use Carp;
use Geo::OGC::Geometry qw/:all/;

our @ISA = qw( Geo::OGC::Geometry );

## @cmethod new(params)
# @brief Construct a new point
# @param params The following syntaxes are allowed:
# @code
# $point = Geo::OGC::Point->new($x, $y);
# $point = Geo::OGC::Point->new($x, $y, $z);
# $point = Geo::OGC::Point->new(point => [$x, $y]);
# $point = Geo::OGC::Point->new(point => [$x, $y, $z]);
# $point = Geo::OGC::Point->new(point => [$x, $y, $z, $m]);
# $point = Geo::OGC::Point->new(pointz => [$x, $y, $z]);
# $point = Geo::OGC::Point->new(pointz => [$x, $y, $z, $m]);
# $point = Geo::OGC::Point->new(pointm => [$x, $y, $m]);
# $point = Geo::OGC::Point->new(pointm => [$x, $y, $z, $m]);
# $point = Geo::OGC::Point->new(pointzm => [$x, $y, $z, $m]);
# $point = Geo::OGC::Point->new(X => $x, Y => $y);
# $point = Geo::OGC::Point->new(X => $x, Y => $y, Z => $z);
# $point = Geo::OGC::Point->new(X => $x, Y => $y, Z => $z, M => $m);
# @endcode
sub new {
    my $package = shift;
    my %params;
    if (@_ == 2 and !($_[0] =~ /^[XYZMpP]/)) { # allow syntax Point->new($x, $y);
	$params{X} = shift;
	$params{Y} = shift;
    } elsif (@_ == 3) { # allow syntax Point->new($x, $y, $z);
	$params{X} = shift;
	$params{Y} = shift;
	$params{Z} = shift;
    } else {
	%params = @_;
    }
    my $self = Geo::OGC::Geometry::new($package, %params);
    return $self;
}

sub init {
    my($self, %params) = @_;
    $self->SUPER::init(%params);
    if ($params{point}) {
	$self->{X} = $params{point}[0];
	$self->{Y} = $params{point}[1];
	$self->{Z} = $params{point}[2] if @{$params{point}} > 2;
	$self->{M} = $params{point}[3] if @{$params{point}} > 3;
    } elsif ($params{pointz}) {
	$self->{X} = $params{pointm}[0];
	$self->{Y} = $params{pointm}[1];
	$self->{Z} = $params{pointm}[2];
	if (@{$params{pointz}} == 4) {
	    $self->{M} = $params{pointm}[3];
	}
    } elsif ($params{pointm}) {
	$self->{X} = $params{pointm}[0];
	$self->{Y} = $params{pointm}[1];
	if (@{$params{pointm}} == 3) {
	    $self->{M} = $params{pointm}[2];
	} elsif (@{$params{pointm}} == 4) {
	    $self->{Z} = $params{pointm}[2];
	    $self->{M} = $params{pointm}[3];
	}
    } elsif ($params{pointzm}) {
	$self->{X} = $params{pointm}[0];
	$self->{Y} = $params{pointm}[1];
	$self->{Z} = $params{pointm}[2];
	$self->{M} = $params{pointm}[3];
    } else {
	$self->{X} = $params{X} if exists $params{X};
	$self->{Y} = $params{Y} if exists $params{Y};
	$self->{Z} = $params{Z} if exists $params{Z};
	$self->{M} = $params{M} if exists $params{M};
    }
}

sub copy {
    my($self, $clone) = @_;
    $self->SUPER::copy($clone);
    for (qw/X Y Z M/) {
	$clone->{$_} = $self->{$_} if exists $self->{$_};
    }
}

## @method point()
# @brief Return a reference to an anonymous array that contains the point data.
# @note Note that there is not difference between [x,y,z] and [x,y,m]
sub point {
    my($self) = @_;
    my @point = ($self->{X}, $self->{Y});
    push @point, $self->{Z} if exists $self->{Z};
    push @point, $self->{M} if exists $self->{M};
    return [@point];
}

sub GeometryType {
    return 'Point';
}

sub Clone {
    my($self) = @_;
    my $m = exists $self->{M} ? 'm' : '';
    return Geo::OGC::Point::new($self, "point$m" => $self->point);
}

sub IsEmpty {
    my($self) = @_;
    return !(exists $self->{X});
}

## @method IsSimple()
# @brief A point is always simple.
sub IsSimple {
    my($self) = @_;
    return 1;
}

sub Is3D {
    my($self) = @_;
    return exists $self->{Z};
}

sub IsMeasured {
    my($self) = @_;
    return exists $self->{M};
}

sub Boundary {
    my($self) = @_;
    return Geo::OGC::GeometryCollection->new();
}

sub X {
    my($self, $X) = @_;
    defined $X ? 
	$self->{X} = $X : $self->{X};
}

sub Y {
    my($self, $Y) = @_;
    defined $Y ? 
	$self->{Y} = $Y : $self->{Y};
}

## @method Z($Z)
# @brief sets or gets the z coordinate
# @param Z [optional]
# @note setting is not in the specification
# @note returns undef if z does not exist or if it exists but is undefined
sub Z {
    my($self, $Z) = @_;
    defined $Z ? 
	$self->{Z} = $Z : (exists $self->{Z} ? $self->{Z} : undef);
}

## @method M($M)
# @brief sets or gets the measure
# @param M [optional]
# @note setting is not in the specification
# @note returns undef if M does not exist or if it exists but is undefined
sub M {
    my($self, $M) = @_;
    defined $M ? 
	$self->{M} = $M : (exists $self->{M} ? $self->{M} : undef);
}

sub as_text {
    my($self, $force_parens, $include_tag) = @_;
    my @coords;
    my $ZM = exists $self->{Z} ? 'Z' : '';
    if (exists $self->{Precision}) {
	for (qw/X Y Z/) {
	    last unless exists $self->{$_};
	    my $s = sprintf('%.'.$self->{Precision}.'e', $self->{$_});
	    push @coords, $s;
	}
    } else {
	for (qw/X Y Z/) {
	    push @coords, $self->{$_} if exists $self->{$_};
	}
    }
    if (exists $self->{M}) {
	push @coords, $self->{M};
	$ZM .= 'M';
    }
    my $text = join(' ', @coords);
    $text = '('.$text.')' if $force_parens;
    $text = "POINT$ZM ".$text if $include_tag;
    return $text;
}

sub Equals {
    my($self, $geom) = @_;
    return 0 unless isa($geom, 'Geo::OGC::Point');
    if (exists $self->{Precision}) {
	for (qw/X Y Z/) {
	    last unless exists $self->{$_} and exists $geom->{$_};
	    my $s = sprintf('%.'.$self->{Precision}.'e', $self->{$_});
	    my $g = sprintf('%.'.$self->{Precision}.'e', $geom->{$_});
	    return 0 if $s != $g;
	}
	return 1;
    }
    # should have an accuracy system that determines the number of meaningful numbers
    return 0 if $self->{X} != $geom->{X};
    return 0 if $self->{Y} != $geom->{Y};
    return 0 if ((exists $self->{Z} and exists $geom->{Z}) and $self->{Z} != $geom->{Z});
    return 1;
}

sub DistanceToLineString {
    my($self, $linestring) = @_;
    my($x, $y) = ($self->{X}, $self->{Y});
    my $distance;
    my $p1;
    for my $p2 (@{$linestring->{Points}}) {
	$p1 = $p2, next unless $p1;
	my $d = distance_point_line($x, $y, $p1->{X}, $p1->{Y}, $p2->{X}, $p2->{Y});
	$distance = $d if !(defined $distance) or $d < $distance;
	$p1 = $p2;
    }
    return $distance;
}

sub Distance {
    my($self, $geom) = @_;
    if (isa($geom, 'Geo::OGC::Point')) {
	return sqrt(($self->{X}-$geom->{X})**2 + ($self->{Y}-$geom->{Y})**2);
    } elsif (isa($geom, 'Geo::OGC::LineString')) {
	return $self->DistanceToLineString($geom);
    } elsif (isa($geom, 'Geo::OGC::Polygon')) {
	if ($geom->{ExteriorRing}->IsPointIn($self)) {
	    for my $ring (@{$geom->{InteriorRings}}) {
		return $ring->DistanceOfPoint($self) if $ring->IsPointIn($self);
	    }
	    return 0;
	} else {
	    return $geom->{ExteriorRing}->DistanceOfPoint($self);
	}
    } elsif (isa($geom, 'Geo::OGC::GeometryCollection')) {
	my $dist = $self->Distance($geom->{Geometries}[0]);
	for my $g (@{$geom->{Geometries}}[1..$#{$geom->{Geometries}}]) {
	    my $d = $self->Distance($g);
	    $dist = $d if $d < $dist;
	}
	return $dist;
    } else {
	croak "can't compute distance between a ".ref($geom)." and a point";
    }
}

# what should this be?
sub Envelope {
    my($self) = @_;
    my $r = Geo::OGC::LinearRing->new;
    $r->AddPoint(Geo::OGC::Point->new($self->{X}, $self->{Y}));
    $r->AddPoint(Geo::OGC::Point->new($self->{X}, $self->{Y}));
    $r->AddPoint(Geo::OGC::Point->new($self->{X}, $self->{Y}));
    $r->AddPoint(Geo::OGC::Point->new($self->{X}, $self->{Y}));
    $r->Close;
    return $r;
}

sub Area {
    return 0;
}

sub MakeCollection {
    my($self) = @_;
    my $coll = Geo::OGC::MultiPoint->new;
    $coll->AddGeometry($self);
    return $coll;
}

sub ApplyTransformation {
    my($self, $transf) = @_;
    if (@_ > 2) {
	($self->{X}, $self->{Y}, $self->{Z}) = $transf->($self->{X}, $self->{Y}, $self->{Z});
    } else {
	($self->{X}, $self->{Y}) = $transf->($self->{X}, $self->{Y});
    }
}

#
#    Curve
#

package Geo::OGC::Curve;
# @brief 1-dimensional geometric object
#
# Curve is implemented as a sequence of Points.

use strict;
use UNIVERSAL qw(isa);
use Carp;
use Geo::OGC::Geometry qw/:all/;

our @ISA = qw( Geo::OGC::Geometry );

sub new {
    my($package, %params) = @_;
    my $self = Geo::OGC::Geometry::new($package, %params);
    return $self;
}

sub init {
    my($self, %params) = @_;
    $self->SUPER::init(%params);
    $self->{Points} = [];
    if ($params{points}) {
	for (@{$params{points}}) {
	    $self->AddPoint(Geo::OGC::Point->new(point=>$_));
	}
    } elsif ($params{pointsm}) {
	for (@{$params{pointsm}}) {
	    $self->AddPoint(Geo::OGC::Point->new(pointm=>$_));
	}
    }
}

sub copy {
    my($self, $clone) = @_;
    $self->SUPER::copy($clone);
    for (@{$self->{Points}}) {
	$clone->AddPoint($_->Clone);
    }
}

sub GeometryType {
    return 'Curve';
}

sub as_text {
    my($self, $force_parens, $include_tag) = @_;
    my $text = join(',', map {$_->as_text} @{$self->{Points}});
    $text = '('.$text.')';
    my $ZM = $self->Is3D ? 'Z' : '';
    $ZM .= 'M' if $self->IsMeasured;
    $text = uc($self->GeometryType).$ZM.' '.$text if $include_tag;
    return $text;
}

## @method AddPoint(point, i)
# @param point A Point object
# @param i [optional] The location in the sequence (1..N+1) where to add the Point. 
# Adds to the end (N+1) by default.
# @note not in the specification
sub AddPoint {
    my($self, $point, $i) = @_;
    croak 'usage: Curve->AddPoint($point) '
	unless $point and isa($point, 'Geo::OGC::Point');
    my $points = $self->{Points};
    if (defined $i) {
	my $temp = $points->[$i-1];
	splice @$points, $i-1, 1, ($point, $temp);
    } else {
	push @$points, $point;
    }
}

## @method DeletePoint(i)
# @param i The location in the sequence (1..N) from where to delete the Point. 
# @note not in the specification
sub DeletePoint {
    my($self, $i) = @_;
    my $points = $self->{Points};
    splice @$points, $i-1, 1;
}

sub StartPoint {
    my($self) = @_;
    my $points = $self->{Points};
    return $points->[0] if @$points;
}

sub EndPoint {
    my($self) = @_;
    my $points = $self->{Points};
    return $points->[$#$points] if @$points;
}

## @method NumPoints()
# @brief Return the number of points in the sequence.
#
# @note Returns all points in a list context.
sub NumPoints {
    my($self) = @_;
    @{$self->{Points}};
}

## @method PointN(N, point)
# @param N A location in the sequence
# @note The first point has the index 1 as OGC SF SQL conformance test uses 1-based indexing. 
# @param point [optional] A Point object, if defined sets the point to index N
sub PointN {
    my($self, $N, $point) = @_;
    my $points = $self->{Points};
    $points->[$N-1] = $point if defined $point;
    return $points->[$N-1];
}

sub Is3D {
    my($self) = @_;
    for (@{$self->{Points}}) {
	return 1 if $_->Is3D;
    }
    return 0;
}

sub IsMeasured {
    my($self) = @_;
    for (@{$self->{Points}}) {
	return 1 if $_->IsMeasured;
    }
    return 0;
}

sub IsClosed {
    my($self) = @_;
    $self->StartPoint()->Equals($self->EndPoint());
}

## @method Close()
# @brief Close the curve by adding the first point also as the last point.
# @note Not in the specification.
sub Close {
    my($self) = @_;
    push @{$self->{Points}}, $self->{Points}[0];
}

## @method IsRing($upgrade)
# @brief Tests whether this curve is a ring, i.e., closed and simple
# @param upgrade [optional, not in the specification] Upgrades this
# curve into a Ring if this really could be a ring
sub IsRing {
    my($self, $upgrade) = @_;
    my $ret = ($self->IsClosed and $self->IsSimple);
    bless($self, 'Geo::OGC::LinearRing') if $ret and $upgrade;
    return $ret;
}

sub Equals {
    my($self, $geom) = @_;
    return 0 unless isa($geom, 'Geo::OGC::Curve');
    return 0 unless $#{$self->{Points}} == $#{$geom->{Points}};
    for my $i (0..$#{$self->{Points}}) {
	return 0 unless $self->{Points}[$i]->Equals($geom->{Points}[$i]);
    }
    return 1;
}

sub Area {
    return 0;
}

sub MakeCollection {
    my($self) = @_;
    my $coll = Geo::OGC::MultiCurve->new;
    $coll->AddGeometry($self);
    return $coll;
}

sub ApplyTransformation {
    my($self, $transf) = @_;
    for my $p (@{$self->{Points}}) {
	$p->ApplyTransformation($transf);
    }
}

## @method Reverse()
# @brief Reverse the order of the points in the sequence.
# @note Not in the specification.
sub Reverse {
    my($self) = @_;
    @{$self->{Points}} = reverse @{$self->{Points}};
}

#
#    LineString
#

package Geo::OGC::LineString;

use strict;
use Carp;
use Geo::OGC::Geometry qw/:all/;

our @ISA = qw( Geo::OGC::Curve );

sub new {
    my($package, %params) = @_;
    my $self = Geo::OGC::Curve::new($package, %params);
    return $self;
}

sub GeometryType {
    return 'LineString';
}

sub IsSimple {
    my($self) = @_;
    my $edges = @{$self->{Points}} - 1;
    return 1 if $edges < 2;
    my $closed = $self->IsClosed;
    my $simple = 1;
    for my $i (0..$edges-1) {
	# check a case of self tangency
	return 0 if $i < $edges-1 and $self->{Points}[$i+2]->Equals($self->{Points}[$i]);
	for my $j ($i+2..$edges-1) {
	    next if $closed and $i == 0 and $j == $edges-1;
	    return 0 if intersect
		(
		 $self->{Points}[$i]{X}, $self->{Points}[$i]{Y},
		 $self->{Points}[$i+1]{X}, $self->{Points}[$i+1]{Y},
		 $self->{Points}[$j]{X}, $self->{Points}[$j]{Y},
		 $self->{Points}[$j+1]{X},$self->{Points}[$j+1]{Y}
		 );
	}
    }
    return 1;
}

sub Envelope {
    my($self) = @_;
    my($minx, $miny, $maxx, $maxy);
    for my $p (@{$self->{Points}}) {
	$minx = $p->{X} if !defined($minx) or $minx > $p->{X};
	$miny = $p->{Y} if !defined($miny) or $miny > $p->{Y};
	$maxx = $p->{X} if !defined($maxx) or $maxx > $p->{X};
	$maxy = $p->{Y} if !defined($maxy) or $maxy > $p->{Y};
    }
    my $r = Geo::OGC::LinearRing->new;
    $r->AddPoint(Geo::OGC::Point->new($minx, $miny));
    $r->AddPoint(Geo::OGC::Point->new($maxx, $miny));
    $r->AddPoint(Geo::OGC::Point->new($maxx, $maxy));
    $r->AddPoint(Geo::OGC::Point->new($minx, $maxy));
    $r->Close;
    return $r;
}

## @method Length()
# @brief The length of this Curve in its associated spatial reference.
# @note Currently computed as a simple euclidean distance.
sub Length {
    my($self) = @_;
    my $l = 0;
    my($x0, $y0) = ($self->{Points}[0]{X}, $self->{Points}[0]{Y});
    for my $i (1..$#{$self->{Points}}) {
	my($x1, $y1) = ($self->{Points}[$i]{X}, $self->{Points}[$i]{Y});
	$l += sqrt(($x1 - $x0)**2+($y1 - $y0)**2);
	($x0, $y0) = ($x1, $y1);
    }
    return $l;
}

sub Distance {
    my($self, $geom) = @_;
    if (isa($geom, 'Geo::OGC::Point')) {
	return $geom->DistanceToLineString($self);
    } elsif (isa($geom, 'Geo::OGC::LineString')) {
	my $dist;
	for my $p (@{$self->{Points}}) {
	    my $d = $p->DistanceToLineString($geom);
	    $dist = $d if !(defined $dist) or $d < $dist;
	}
	return $dist;
    } elsif (isa($geom, 'Geo::OGC::Polygon')) {
	my $dist;
	for my $p (@{$self->{Points}}) {
	    my $d = $p->Distance($geom);
	    $dist = $d if !(defined $dist) or $d < $dist;
	}
	return $dist;
    } elsif (isa($geom, 'Geo::OGC::GeometryCollection')) {
	my $dist = $self->Distance($geom->{Geometries}[0]);
	for my $g (@{$geom->{Geometries}}[1..$#{$geom->{Geometries}}]) {
	    my $d = $self->Distance($g);
	    $dist = $d if $d < $dist;
	}
	return $dist;
    } else {
	croak "can't compute distance between a ".ref($geom)." and a line string";
    }
}

sub MakeCollection {
    my($self) = @_;
    my $coll = Geo::OGC::MultiLineString->new;
    $coll->AddGeometry($self);
    return $coll;
}

#
#    Line
#

package Geo::OGC::Line;

use strict;
use Carp;

our @ISA = qw( Geo::OGC::LineString );

sub new {
    my($package, %params) = @_;
    my $self = Geo::OGC::LineString::new($package, %params);
    return $self;
}

sub GeometryType {
    return 'Line';
}

#
#    LinearRing
#

package Geo::OGC::LinearRing;

use strict;
use Carp;
use Geo::OGC::Geometry qw/:all/;

our @ISA = qw( Geo::OGC::LineString );

sub new {
    my($package, %params) = @_;
    my $self = Geo::OGC::LineString::new($package, %params);
    return $self;
}

sub GeometryType {
    return 'LinearRing';
}

## @method IsPointIn(Geo::OGC::Point point)
# @brief Tests whether the given point is within the ring
# @note uses the pnpoly algorithm from
# http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
# @note Assumes a simple closed ring
sub IsPointIn {
    my($self, $point) = @_;
    my($x, $y) = ($point->{X}, $point->{Y});
    my $c = 0;
    my $prev;
    for (@{$self->{Points}}) {
	$prev = $_, next unless $prev;
        $c = !$c if (((( $_->{Y} <= $y ) && ( $y < $prev->{Y} )) ||
		      (( $prev->{Y} <= $y ) && ( $y < $_->{Y} ))) &&
		     ( $x < ( $prev->{X} - $_->{X} ) * 
		       ( $y - $_->{Y} ) / ( $prev->{Y} - $_->{Y} ) + $_->{X} ));
	$prev = $_;
    }
    return $c;
}

## @method Area()
# @brief Computes the area of the ring.
# @note Not in the specification.
# @note Computed as a simple euclidean area.
# @note Assumes a simple closed ring
# @return area or negative area if the sense of the rotation of the ring is clockwise
sub Area {
    my($self) = @_;
    my $area = 0;
    my $N = $#{$self->{Points}}-1; # skip the closing point
    my $j = 0;
    for my $i (0..$N) {
	$j++;
	$area += $self->{Points}[$i]{X} * $self->{Points}[$j]{Y};
	$area -= $self->{Points}[$i]{Y} * $self->{Points}[$j]{X};
    }
    return $area/2;
}

sub Centroid {
    my($self) = @_;
    my($area, $x, $y) = (0, 0, 0);
    my $N = $#{$self->{Points}}-1; # skip the closing point
    for my $i (0..$N) {
	my($xi, $yi, $xj, $yj) = ( $self->{Points}[$i]{X}, $self->{Points}[$i]{Y},
				   $self->{Points}[$i+1]{X}, $self->{Points}[$i+1]{Y} );
	my $b = $xi * $yj - $xj * $yi;
	$area += $b;
	$x += ($xi + $xj) * $b;
	$y += ($yi + $yj) * $b;
    }
    $x /= abs(3*$area); # 6 but $area is 2 * real area
    $y /= abs(3*$area);
    return Geo::OGC::Point->new($x, $y);
}

## @method IsCCW()
# @brief Returns true if the points in this ring are arranged counterclockwise.
# @note Not in the specification.
# @note Assumes a simple closed ring.
sub IsCCW {
    my($self) = @_;
    # find the northernmost point
    my $t = 0;
    my $N = $#{$self->{Points}}-1; # skip the closing point
    for my $i (1..$N) {
	$t = $i if $self->{Points}[$i]{Y} > $self->{Points}[$t]{Y};
    }
    # the previous point
    my $p = $t-1;
    $p = $N if $p < 0;
    # the next point
    my $n = $t+1;
    $n = 0 if $n > $N;
    return ccw($self->{Points}[$p]{X}, $self->{Points}[$p]{Y}, 
	       $self->{Points}[$t]{X}, $self->{Points}[$t]{Y}, 
	       $self->{Points}[$n]{X}, $self->{Points}[$n]{Y}) == 1;
}

#
#    Surface
#

package Geo::OGC::Surface;

use strict;
use Carp;

our @ISA = qw( Geo::OGC::Geometry );

sub new {
    my($package, %params) = @_;
    my $self = Geo::OGC::Geometry::new($package, %params);
    return $self;
}

sub GeometryType {
    return 'Surface';
}

sub Area {
    my($self) = @_;
    croak "Area method for class ".ref($self)." is not implemented yet";
}

sub Centroid {
    my($self) = @_;
    croak "Centroid method for class ".ref($self)." is not implemented yet";
}

sub PointOnSurface {
    my($self) = @_;
    croak "PointOnSurface method for class ".ref($self)." is not implemented yet";
}

sub MakeCollection {
    my($self) = @_;
    my $coll = Geo::OGC::MultiSurface->new;
    $coll->AddGeometry($self);
    return $coll;
}

#
#    Polygon
#

package Geo::OGC::Polygon;

use strict;
use UNIVERSAL qw(isa);
use Carp;

our @ISA = qw( Geo::OGC::Surface );

sub new {
    my($package, %params) = @_;
    my $self = Geo::OGC::Surface::new($package, %params);
    return $self;
}

sub init {
    my($self, %params) = @_;
    $self->SUPER::init(%params);
    $self->{ExteriorRing} = undef;
    $self->{InteriorRings} = [];
}

sub copy {
    my($self, $clone) = @_;
    $self->SUPER::copy($clone);
    $clone->ExteriorRing($self->{ExteriorRing}->Clone) if $self->{ExteriorRing};
    for (@{$self->{InteriorRings}}) {
	$clone->AddInteriorRing($_->Clone);
    }
}

sub GeometryType {
    return 'Polygon';
}

sub Is3D {
    my($self) = @_;
    return 1 if $self->{ExteriorRing}->Is3D;
    for (@{$self->{InteriorRings}}) {
	return 1 if $_->Is3D;
    }
    return 0;
}

sub IsMeasured {
    my($self) = @_;
    return 1 if $self->{ExteriorRing}->IsMeasured;
    for (@{$self->{InteriorRings}}) {
	return 1 if $_->IsMeasured;
    }
    return 0;
}

sub AddInteriorRing {
    my($self, $ring, $i) = @_;
    croak 'usage: Polygon->AddInteriorRing($ring[, $i])' 
	unless $ring and isa($ring, 'Geo::OGC::LinearRing');
    my $rings = $self->{InteriorRings};
    $i = @$rings unless defined $i;
    if (@$rings) {
	my $temp = $rings->[$i-1];
	splice @$rings,$i-1,1,($temp, $ring);
    } else {
	push @$rings, $ring;
    }
}

sub ExteriorRing {
    my($self, $ring) = @_;
    if (defined $ring) {
	croak 'usage: Polygon->ExteriorRing($ring)' 
	    unless isa($ring, 'Geo::OGC::LinearRing');
	$self->{ExteriorRing} = $ring;
    } else {
	return $self->{ExteriorRing};
    }
}

sub Envelope {
    my($self) = @_;
    return $self->{ExteriorRing}->Envelope;
}

## @method NumInteriorRing()
# @brief Return the number of interior rings in the polygon.
#
# @note Returns all interior rings in a list context.
sub NumInteriorRing {
    my($self) = @_;
    @{$self->{InteriorRings}};
}

sub InteriorRingN {
    my($self, $N, $ring) = @_;
    my $rings = $self->{InteriorRings};
    $rings->[$N-1] = $ring if defined $ring;
    return $rings->[$N-1] if @$rings;
}

# @note Assumes the order of the interior rings is the same.
sub Equals {
    my($self, $geom) = @_;
    return 0 unless isa($geom, 'Geo::OGC::Polygon');
    return 0 unless @{$self->{InteriorRings}} == @{$geom->{InteriorRings}};
    return 0 unless $self->{ExteriorRing}->Equals($geom->{ExteriorRing});
    for my $i (0..$#{$self->{InteriorRings}}) {
	return 0 unless $self->{InteriorRings}[$i]->Equals($geom->{InteriorRings}[$i]);
    }
    return 1;
}

sub Area {
    my($self) = @_;
    my $a = $self->{ExteriorRing}->Area;
    for my $ring (@{$self->{InteriorRings}}) {
	$a -= $ring->Area;
    }
    return $a;
}

sub IsPointIn {
    my($self, $point) = @_;
    my $c = $self->{ExteriorRing}->IsPointIn($point);
    if ($c) {
	for my $ring (@{$self->{InteriorRings}}) {
	    $c = 0, last if $ring->IsPointIn($point);
	}
    }
    return $c;
}

sub as_text {
    my($self, $force_parens, $include_tag) = @_;
    my $text .= $self->{ExteriorRing}->as_text if $self->{ExteriorRing};
    for my $ring (@{$self->{InteriorRings}}) {
	$text .= ', ';
	$text .= $ring->as_text;
    }
    my $ZM = $self->Is3D ? 'Z' : '';
    $ZM .= 'M' if $self->IsMeasured;
    $text = uc($self->GeometryType).$ZM.' ('.$text.')' if $include_tag;
    return $text;
}

sub MakeCollection {
    my($self) = @_;
    my $coll = Geo::OGC::MultiPolygon->new;
    $coll->AddGeometry($self);
    return $coll;
}

#
#    Triangle
#

package Geo::OGC::Triangle;

use strict;
use UNIVERSAL qw(isa);
use Carp;

our @ISA = qw( Geo::OGC::Polygon );

sub new {
    my($package, %params) = @_;
    my $self = Geo::OGC::Polygon::new($package, %params);
    return $self;
}

sub GeometryType {
    return 'Triangle';
}

#
#    PolyhedralSurface
#

package Geo::OGC::PolyhedralSurface;

use strict;
use UNIVERSAL qw(isa);
use Carp;

our @ISA = qw( Geo::OGC::Surface );

sub new {
    my($package, %params) = @_;
    my $self = Geo::OGC::Surface::new($package, %params);
    return $self;
}

sub init {
    my($self, %params) = @_;
    $self->SUPER::init(%params);
    $self->{Patches} = []; # linearrings
    if ($params{patches}) {
	for (@{$params{patches}}) {
	    $self->AddPatch(Geo::OGC::LinearRing->new(points=>$_));
	}
    } elsif ($params{patchesm}) {
	for (@{$params{patches}}) {
	    $self->AddPatch(Geo::OGC::LinearRing->new(pointsm=>$_));
	}
    }
}

sub copy {
    my($self, $clone) = @_;
    $self->SUPER::copy($clone);
    for (@{$self->{Patches}}){
	$clone->AddPatch($_->Clone);
    }
}

sub GeometryType {
    return 'PolyhedralSurface';
}

sub AddPatch {
    my($self, $patch, $i) = @_;
    croak 'usage: PolyhedralSurface->AddPatch($patch[, $i])' 
	unless $patch and isa($patch, 'Geo::OGC::LinearRing');
    my $patches = $self->{Patches};
    $i = @$patches unless defined $i;
    if (@$patches) {
	my $temp = $patches->[$i-1];
	splice @$patches,$i-1,1,($temp, $patch);
    } else {
	push @$patches, $patch;
    }
}

sub NumPatches {
    my($self) = @_;
    @{$self->{Patches}};
}

sub PatchN {
    my($self, $N, $patch) = @_;
    my $patches = $self->{Patches};
    $patches->[$N-1] = $patch if defined $patch;
    return $patches->[$N-1] if @$patches;
}

sub BoundingPolygons {
    my($self, $p) = @_;
    croak "BoundingPolygons method for class ".ref($self)." is not implemented yet";
}

sub IsClosed {
    my($self) = @_;
    croak "IsClosed method for class ".ref($self)." is not implemented yet";
}

sub IsMeasured {
    my($self) = @_;
    for (@{$self->{Patches}}) {
	return 1 if $_->IsMeasured;
    }
    return 0;
}

sub as_text {
    my($self, $force_parens, $include_tag) = @_;
    my $text = '';
    for my $patch (@{$self->{Patches}}) {
	$text .= ', ';
	$text .= $patch->as_text;
    }
    my $ZM = $self->Is3D ? 'Z' : '';
    $ZM .= 'M' if $self->IsMeasured;
    $text = uc($self->GeometryType).$ZM.' ('.$text.')' if $include_tag;
    return $text;
}

#
#    TIN
#

package Geo::OGC::TIN;

use strict;
use UNIVERSAL qw(isa);
use Carp;

our @ISA = qw( Geo::OGC::PolyhedralSurface );

sub new {
    my($package, %params) = @_;
    my $self = Geo::OGC::Surface::new($package, %params);
    return $self;
}

sub GeometryType {
    return 'TIN';
}

#
#    GeometryCollection
#

package Geo::OGC::GeometryCollection;

use strict;
use UNIVERSAL qw(isa);
use Carp;

our @ISA = qw( Geo::OGC::Geometry );

sub new {
    my($package, %params) = @_;
    my $self = Geo::OGC::Geometry::new($package, %params);
    return $self;
}

sub init {
    my($self, %params) = @_;
    $self->SUPER::init(%params);
    $self->{Geometries} = [];
}

sub copy {
    my($self, $clone) = @_;
    $self->SUPER::copy($clone);
    for (@{$self->{Geometries}}) {
	$clone->AddGeometry($_->Clone);
    }
}

sub GeometryType {
    return 'GeometryCollection';
}

sub Is3D {
    my($self) = @_;
    for (@{$self->{Geometries}}) {
	return 1 if $_->Is3D;
    }
    return 0;
}

sub IsMeasured {
    my($self) = @_;
    for (@{$self->{Geometries}}) {
	return 1 if $_->IsMeasured;
    }
    return 0;
}

sub as_text {
    my($self, $force_parens, $include_tag) = @_;
    my $text = join(',', map {$_->as_text(1, 1)} @{$self->{Geometries}});
    my $ZM = $self->Is3D ? 'Z' : '';
    $ZM .= 'M' if $self->IsMeasured;
    $text = uc($self->GeometryType).$ZM.' ('.$text.')' if $include_tag;
    return $text;
}

sub ElementType {
    return 'Geometry';
}

sub AddGeometry {
    my($self, $geometry, $i) = @_;
    croak 'usage: GeometryCollection->AddGeometry($geometry[, $i])' 
	unless $geometry and isa($geometry, 'Geo::OGC::Geometry');
    my $geometries = $self->{Geometries};
    $i = @$geometries unless defined $i;
    if (@$geometries) {
	my $temp = $geometries->[$i-1];
	splice @$geometries,$i-1,1,($temp, $geometry);
    } else {
	push @$geometries, $geometry;
    }
}

## @method NumGeometries()
# @brief Return the number of geometries in the collection.
#
# @note Returns all geometries in a list context.
sub NumGeometries {
    my($self) = @_;
    @{$self->{Geometries}};
}

## @method GeometryN(N)
# @brief Return the Nth geometry from the collection 
# (the index of the first geometry is 1).
#
sub GeometryN {
    my($self, $N, $geometry) = @_;
    my $geometries = $self->{Geometries};
    $geometries->[$N-1] = $geometry if defined $geometry;
    return $geometries->[$N-1] if @$geometries;
}

sub Envelope {
    my($self) = @_;
    my($minx, $miny, $maxx, $maxy);
    for my $g (@{$self->{Geometries}}) {
	my $e = $g->Envelope;
	my $min = $e->PointN(1);
	my $max = $e->PointN(3);
	$minx = $min->{X} if !defined($minx) or $minx > $min->{X};
	$miny = $min->{Y} if !defined($miny) or $miny > $min->{Y};
	$maxx = $max->{X} if !defined($maxx) or $maxx > $max->{X};
	$maxy = $max->{Y} if !defined($maxy) or $maxy > $max->{Y};
    }
    my $r = new Geo::OGC::LinearRing;
    $r->AddPoint(Geo::OGC::Point->new($minx, $miny));
    $r->AddPoint(Geo::OGC::Point->new($maxx, $miny));
    $r->AddPoint(Geo::OGC::Point->new($maxx, $maxy));
    $r->AddPoint(Geo::OGC::Point->new($minx, $maxy));
    $r->Close;
    return $r;
}

# @note Assumes the order is the same.
sub Equals {
    my($self, $geom) = @_;
    return 0 unless isa($geom, 'Geo::OGC::GeometryCollection');
    return 0 unless @{$self->{Geometries}} == @{$geom->{Geometries}};
    for my $i (0..$#{$self->{Geometries}}) {
	return 0 unless $self->{Geometries}[$i]->Equals($geom->{Geometries}[$i]);
    }
    return 1;
}

#
#    MultiSurface
#

package Geo::OGC::MultiSurface;

use strict;
use Carp;

our @ISA = qw( Geo::OGC::GeometryCollection );

sub new {
    my($package, %params) = @_;
    my $self = Geo::OGC::GeometryCollection::new($package, %params);
    return $self;
}

sub GeometryType {
    return 'MultiSurface';
}

sub ElementType {
    return 'Surface';
}

sub Area {
    my($self) = @_;
    croak "Area method for class ".ref($self)." is not implemented yet";
}

sub Centroid {
    my($self) = @_;
    croak "Centroid method for class ".ref($self)." is not implemented yet";
}

sub PointOnSurface {
    my($self) = @_;
    croak "PointOnSurface method for class ".ref($self)." is not implemented yet";
}

#
#    MultiCurve
#

package Geo::OGC::MultiCurve;

use strict;
use Carp;

our @ISA = qw( Geo::OGC::GeometryCollection );

sub new {
    my($package, %params) = @_;
    my $self = Geo::OGC::GeometryCollection::new($package, %params);
    return $self;
}

sub GeometryType {
    return 'MultiCurve';
}

sub ElementType {
    return 'Curve';
}

#
#    MultiPoint
#

package Geo::OGC::MultiPoint;

use strict;
use Carp;

our @ISA = qw( Geo::OGC::GeometryCollection );

sub new {
    my($package, %params) = @_;
    my $self = Geo::OGC::GeometryCollection::new($package, %params);
    return $self;
}

sub init {
    my($self, %params) = @_;
    $self->SUPER::init(%params);
    if ($params{points}) {
	for (@{$params{points}}) {
	    $self->AddGeometry(Geo::OGC::Point->new(point=>$_));
	}
    } elsif ($params{pointsm}) {
	for (@{$params{pointsm}}) {
	    $self->AddGeometry(Geo::OGC::Point->new(pointm=>$_));
	}
    }
}

sub GeometryType {
    return 'MultiPoint';
}

sub as_text {
    my($self, $force_parens, $include_tag) = @_;
    my $text = join(',', map {$_->as_text(1)} @{$self->{Geometries}});
    my $ZM = $self->Is3D ? 'Z' : '';
    $ZM .= 'M' if $self->IsMeasured;
    $text = uc($self->GeometryType).$ZM.' ('.$text.')' if $include_tag;
    return $text;
}

sub ElementType {
    return 'Point';
}

sub Boundary {
    my($self) = @_;
    return Geo::OGC::GeometryCollection->new();
}

#
#    MultiPolygon
#

package Geo::OGC::MultiPolygon;

use strict;
use Carp;

our @ISA = qw( Geo::OGC::MultiSurface );

sub new {
    my($package, %params) = @_;
    my $self = Geo::OGC::MultiSurface::new($package, %params);
    return $self;
}

sub GeometryType {
    return 'MultiPolygon';
}

sub as_text {
    my($self, $force_parens, $include_tag) = @_;
    my $text = join(',', map {$_->as_text(1)} @{$self->{Geometries}});
    my $ZM = $self->Is3D ? 'Z' : '';
    $ZM .= 'M' if $self->IsMeasured;
    $text = uc($self->GeometryType).$ZM.' ('.$text.')' if $include_tag;
    return $text;
}

sub ElementType {
    return 'Polygon';
}

#
#    MultiLineString
#

package Geo::OGC::MultiLineString;

use strict;
use Carp;

our @ISA = qw( Geo::OGC::MultiCurve );

sub new {
    my($package, %params) = @_;
    my $self = Geo::OGC::MultiCurve::new($package, %params);
    return $self;
}

sub GeometryType {
    return 'MultiLineString';
}

sub as_text {
    my($self, $force_parens, $include_tag) = @_;
    my $text = join(',', map {$_->as_text(1)} @{$self->{Geometries}});
    my $ZM = $self->Is3D ? 'Z' : '';
    $ZM .= 'M' if $self->IsMeasured;
    $text = uc($self->GeometryType).$ZM.' ('.$text.')' if $include_tag;
    return $text;
}

sub ElementType {
    return 'LineString';
}
