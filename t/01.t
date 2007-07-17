# Methods of Curve

use UNIVERSAL qw(isa);
use Test::More qw(no_plan);

BEGIN { 
    use_ok( 'Geo::OGC::Geometry' );
}

$c = new Geo::OGC::Curve;
$c->AddPoint(Geo::OGC::Point->new(0, 0));
$c->AddPoint(Geo::OGC::Point->new(1, 0));
$c->AddPoint(Geo::OGC::Point->new(1, 1));
$c->AddPoint(Geo::OGC::Point->new(0, 1));

$c->AddPoint(Geo::OGC::Point->new(2, 0), 1);
ok($c->PointN(1)->Equals(Geo::OGC::Point->new(2, 0)), "AddPoint 1 a");
ok($c->PointN(2)->Equals(Geo::OGC::Point->new(0, 0)), "AddPoint 1 b");
$c->DeletePoint(1);
ok($c->PointN(1)->Equals(Geo::OGC::Point->new(0, 0)), "DeletePoint 1");

$c->AddPoint(Geo::OGC::Point->new(2, 0), 2);
ok($c->PointN(2)->Equals(Geo::OGC::Point->new(2, 0)), "AddPoint 2 a");
ok($c->PointN(3)->Equals(Geo::OGC::Point->new(1, 0)), "AddPoint 2 b");
$c->DeletePoint(2);
ok($c->PointN(2)->Equals(Geo::OGC::Point->new(1, 0)), "DeletePoint 2");

$c->AddPoint(Geo::OGC::Point->new(2, 0), 4);
ok($c->PointN(4)->Equals(Geo::OGC::Point->new(2, 0)), "AddPoint 4 a");
ok($c->PointN(5)->Equals(Geo::OGC::Point->new(0, 1)), "AddPoint 4 b");
$c->DeletePoint(4);
ok($c->PointN(4)->Equals(Geo::OGC::Point->new(0, 1)), "DeletePoint 4");

bless $c, 'Geo::OGC::LineString';
ok($c->Length == 3, "Length");
ok($c->StartPoint->Equals(Geo::OGC::Point->new(0, 0)), "StartPoint");
ok($c->EndPoint->Equals(Geo::OGC::Point->new(0, 1)), "EndPoint");
ok($c->NumPoints == 4, "NumPoints");
ok(!$c->Is3D, "Is3D a");
$c->AddPoint(Geo::OGC::Point->new(2, 0, 3), 1);
ok($c->Is3D, "Is3D b");
$c->DeletePoint(1);

ok($c->IsSimple, "IsSimple a");
$c->AddPoint(Geo::OGC::Point->new(0.5, -1));
ok(!$c->IsSimple, "IsSimple b");
$c->DeletePoint(5);
ok($c->IsSimple, "IsSimple c");

ok(!$c->IsClosed, "IsClosed a");
$c->Close;
ok($c->IsClosed, "IsClosed b");

ok($c->IsRing(1), "IsRing");
ok($c->Area == 1, "Area");

$c = new Geo::OGC::LineString;
$c->AddPoint(Geo::OGC::Point->new(0, 0));
$c->AddPoint(Geo::OGC::Point->new(1, 0));
$c->AddPoint(Geo::OGC::Point->new(1, 1));
$c->Close;
ok($c->IsClosed, "IsClosed 3");
ok($c->IsSimple, "IsSimple 3");
ok($c->IsRing(1), "IsRing 3");
my @a = $c->Area;
ok($a[0] == 0.5, "Area 3");
$c->Reverse;
@a = $c->Area;
ok($a[0] == -0.5, "Area 3 b");

