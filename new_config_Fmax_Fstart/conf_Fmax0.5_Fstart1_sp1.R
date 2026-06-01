# Eco-OSMOSE-NS MAIN CONFIGURATION FILE
# North Sea
# OSMOSE version 4.4.0

# All parameters ---------------------------------------------------------------
osmose.version = 4.4.1

mortality.algorithm = stochastic
mortality.subdt = 10
# stochastic.mortality.seed = 21
# movement.randomseed.fixed = FALSE

## Species list:

# Focus:
#   0-	Atlantic Herring (Clupea harengus)
#   1-	Atlantic Mackerel (Scomber scombrus)
#   2-	Lesser sand-eel (Ammodytes marinus)
#   3-	European sprat (Sprattus sprattus)
#   4-  Norway pout (Trisopterus esmarkii)
#   5-  European plaice (Pleuronectes platessa)
#   6-  Common sole (Solea solea)
#   7-  Saithe (Pollachius virens)
#   8-  Atlantic cod (Gadus morhua)
#   9-  Haddock (Melanogrammus aeglefinus)
#   10-	Atlantic horse mackerel (Trachurus trachurus)
#   11-	Whiting (Merlangius merlangus)
#   12- Common dab (Limanda limanda)
#   13- Grey gurnard (Eutrigla gurnardus)
#   14-	European hake (Merluccius merluccius)
#   15-	Shrimp (Crangon crangon)
# Resource:
#   16-	Microphytoplankton
#   17-	Diatoms
#   18-	Heterotrophic flagellates
#   19-	Microzooplankton
#   20-	Mesozooplankton
#   21- Suspension feeders
#   22- Deposit feeders
#   23- Meiobenthos
#   24- Large benthos
#   25- Very large benthos

# Calibration parameters -------------------------------------------------------

# mortality.additional.larva.rate.sp0 = 6.10349484697158
# mortality.additional.larva.rate.sp1 = 4.64743912083943
# mortality.additional.larva.rate.sp2 = 2.25772421006519
# mortality.additional.larva.rate.sp3 = 0.93653180427938
# mortality.additional.larva.rate.sp4 = 3.76131651995884
# mortality.additional.larva.rate.sp5 = 5.21619885928796
# mortality.additional.larva.rate.sp6 = 7.01498455642933
# mortality.additional.larva.rate.sp7 = 2.72955476392776
# mortality.additional.larva.rate.sp8 = 9.34281126685825
# mortality.additional.larva.rate.sp9 = 2.82198045352037
# mortality.additional.larva.rate.sp10 = 1.20470557168576
# mortality.additional.larva.rate.sp11 = 8.56133204806727
# mortality.additional.larva.rate.sp12 = 4.63271893188496
# mortality.additional.larva.rate.sp13 = 5.15337726111462
# mortality.additional.larva.rate.sp14 = 7.7196438205665
# mortality.additional.larva.rate.sp15 = 0.0572525600128366

# new larval mortality (x ndt, ndt=24)
mortality.additional.larva.rate.sp0 = 146.483876327318
mortality.additional.larva.rate.sp1 = 111.538538900146
mortality.additional.larva.rate.sp2 = 54.1853810415646
mortality.additional.larva.rate.sp3 = 22.4767633027051
mortality.additional.larva.rate.sp4 = 90.2715964790122
mortality.additional.larva.rate.sp5 = 125.188772622911
mortality.additional.larva.rate.sp6 = 168.359629354304
mortality.additional.larva.rate.sp7 = 65.5093143342662
mortality.additional.larva.rate.sp8 = 224.227470404598
mortality.additional.larva.rate.sp9 = 67.7275308844889
mortality.additional.larva.rate.sp10 = 28.9129337204582
mortality.additional.larva.rate.sp11 = 205.471969153614
mortality.additional.larva.rate.sp12 = 111.185254365239
mortality.additional.larva.rate.sp13 = 123.681054266751
mortality.additional.larva.rate.sp14 = 185.271451693596
mortality.additional.larva.rate.sp15 = 1.37406144030808

species.accessibility2fish.sp16 = 0.132785620764358
species.accessibility2fish.sp17 = 0.0281234083561995
species.accessibility2fish.sp18 = 0.186714399330888
species.accessibility2fish.sp19 = 0.0321694042977577
species.accessibility2fish.sp20 = 0.0738045513264178
species.accessibility2fish.sp21 = 0.00183914566903387
species.accessibility2fish.sp22 = 5.71676156957018e-05
species.accessibility2fish.sp23 = 0.000939821640504932
species.accessibility2fish.sp24 = 0.00939960490393895
species.accessibility2fish.sp25 = 0.0145799271877107

predation.ingestion.rate.max.sp0 = 14.0161771995127
predation.ingestion.rate.max.sp1 = 16.4430568215567
predation.ingestion.rate.max.sp2 = 9.48185907748328
predation.ingestion.rate.max.sp3 = 12.119810661574
predation.ingestion.rate.max.sp4 = 9.91981236341634
predation.ingestion.rate.max.sp5 = 10.9236872467061
predation.ingestion.rate.max.sp6 = 9.94816846789243
predation.ingestion.rate.max.sp7 = 14.4507137847945
predation.ingestion.rate.max.sp8 = 21.5633912243199
predation.ingestion.rate.max.sp9 = 15.7906712579798
predation.ingestion.rate.max.sp10 = 13.3885972672411
predation.ingestion.rate.max.sp11 = 17.8426721069184
predation.ingestion.rate.max.sp12 = 8.78878491915743
predation.ingestion.rate.max.sp13 = 13.5696659330434
predation.ingestion.rate.max.sp14 = 17.2822647633446
predation.ingestion.rate.max.sp15 = 8.91297095163207

F.factor.sp0 = 1.26580601110599
F.factor.sp1 = 2.56100968753689
F.factor.sp2 = 1.61952301469531
F.factor.sp3 = 1.33095773934421
F.factor.sp4 = 1.38839713307704
F.factor.sp5 = 1.14675509827372
F.factor.sp6 = 0.96115620294063
F.factor.sp7 = 0.968361758445048
F.factor.sp8 = 1.12545641380261
F.factor.sp9 = 1.36592394280662
F.factor.sp10 = 1.26545734889451
F.factor.sp11 = 1.86283048571559
F.factor.sp12 = 1.63107970606178
F.factor.sp13 = 1.62020777222316
F.factor.sp14 = 1.22196472076686
F.factor.sp15 = 0.425135134132864

mortality.additional.rate.sp0 = 0.11601151182773
mortality.additional.rate.sp1 = 0.251529724669577
mortality.additional.rate.sp2 = 0.138320973545911
mortality.additional.rate.sp3 = 0.0702435464333923
mortality.additional.rate.sp4 = 0.263994247897901
mortality.additional.rate.sp5 = 0.0177832969813442
mortality.additional.rate.sp6 = 0.298005973278076
mortality.additional.rate.sp7 = 0.550526966152144
mortality.additional.rate.sp8 = 0.527649887669659
mortality.additional.rate.sp9 = 0.75954600216465
mortality.additional.rate.sp10 = 0.2299911889241
mortality.additional.rate.sp11 = 0.10992469807082
mortality.additional.rate.sp12 = 0.200789223641564
mortality.additional.rate.sp13 = 0.0265415982712645
mortality.additional.rate.sp14 = 0.271740155170216
mortality.additional.rate.sp15 = 0.160236121241078

# Grid parameters --------------------------------------------------------------

grid.netcdf.file = LTL_ns_corrected_one.nc
grid.var.lat = ny
grid.var.lon = nx
grid.var.mask = Diatoms

# Species parameters -----------------------------------------------------------

species.egg.size.sp0 = 0
species.egg.size.sp1 = 0
species.egg.size.sp2 = 10
species.egg.size.sp3 = 0
species.egg.size.sp4 = 0
species.egg.size.sp5 = 0
species.egg.size.sp6 = 0
species.egg.size.sp7 = 0
species.egg.size.sp8 = 0
species.egg.size.sp9 = 0
species.egg.size.sp10 = 0
species.egg.size.sp11 = 0
species.egg.size.sp12 = 0
species.egg.size.sp13 = 0
species.egg.size.sp14 = 0
species.egg.size.sp15 = 0

species.egg.weight.sp0 = 0.000667459
species.egg.weight.sp1 = 0.000585947
species.egg.weight.sp2 = 0.000997299
species.egg.weight.sp3 = 0.000280466
species.egg.weight.sp4 = 0.000157231
species.egg.weight.sp5 = 0.002038534
species.egg.weight.sp6 = 0.000581087
species.egg.weight.sp7 = 0.000330484
species.egg.weight.sp8 = 0.000455192
species.egg.weight.sp9 = 0.001361162
species.egg.weight.sp10 = 0.000735224
species.egg.weight.sp11 = 0.000428647
species.egg.weight.sp12 = 0.000172213
species.egg.weight.sp13 = 0.000992277
species.egg.weight.sp14 = 0.000672383
species.egg.weight.sp15 = 0.002031007

species.lifespan.sp0 = 17
species.lifespan.sp1 = 17
species.lifespan.sp2 = 10
species.lifespan.sp3 = 10
species.lifespan.sp4 = 7
species.lifespan.sp5 = 25
species.lifespan.sp6 = 17
species.lifespan.sp7 = 21
species.lifespan.sp8 = 17
species.lifespan.sp9 = 18
species.lifespan.sp10 = 15
species.lifespan.sp11 = 16
species.lifespan.sp12 = 10
species.lifespan.sp13 = 15
species.lifespan.sp14 = 15
species.lifespan.sp15 = 3

species.name.sp0 = AtlanticHerring
species.name.sp1 = AtlanticMackerel
species.name.sp2 = Sandeel
species.name.sp3 = EuropeanSprat
species.name.sp4 = NorwayPout
species.name.sp5 = EuropeanPlaice
species.name.sp6 = CommonSole
species.name.sp7 = Saithe
species.name.sp8 = AtlanticCod
species.name.sp9 = Haddock
species.name.sp10 = HorseMackerel
species.name.sp11 = Whiting
species.name.sp12 = CommonDab
species.name.sp13 = GreyGurnard
species.name.sp14 = Hake
species.name.sp15 = Shrimp

species.relativefecundity.sp0 = 1159
species.relativefecundity.sp1 = 1275
species.relativefecundity.sp2 = 1202
species.relativefecundity.sp3 = 1993
species.relativefecundity.sp4 = 3199
species.relativefecundity.sp5 = 255
species.relativefecundity.sp6 = 1040
species.relativefecundity.sp7 = 630
species.relativefecundity.sp8 = 1188
species.relativefecundity.sp9 = 430
species.relativefecundity.sp10 = 1655
species.relativefecundity.sp11 = 1700
species.relativefecundity.sp12 = 2931
species.relativefecundity.sp13 = 873
species.relativefecundity.sp14 = 426
species.relativefecundity.sp15 = 815

species.sexratio.sp0 = 0.5
species.sexratio.sp1 = 0.5
species.sexratio.sp2 = 0.5
species.sexratio.sp3 = 0.5
species.sexratio.sp4 = 0.5
species.sexratio.sp5 = 0.5
species.sexratio.sp6 = 0.5
species.sexratio.sp7 = 0.5
species.sexratio.sp8 = 0.5
species.sexratio.sp9 = 0.5
species.sexratio.sp10 = 0.5
species.sexratio.sp11 = 0.5
species.sexratio.sp12 = 0.5
species.sexratio.sp13 = 0.5
species.sexratio.sp14 = 0.5
species.sexratio.sp15 = 0.5

species.K.sp0 = 0.44
species.K.sp1 = 0.837
species.K.sp2 = 0.36
species.K.sp3 = 0.503
species.K.sp4 = 0.52
species.K.sp5 = 0.16
species.K.sp6 = 0.148
species.K.sp7 = 0.211
species.K.sp8 = 0.208
species.K.sp9 = 0.25
species.K.sp10 = 0.24
species.K.sp11 = 0.432
species.K.sp12 = 0.21
species.K.sp13 = 0.16
species.K.sp14 = 0.073
species.K.sp15 = 0.208

species.length2weight.allometric.power.sp0 = 3.40588401
species.length2weight.allometric.power.sp1 = 3.062840519
species.length2weight.allometric.power.sp2 = 3.44
species.length2weight.allometric.power.sp3 = 3.148474632
species.length2weight.allometric.power.sp4 = 3.096696558
species.length2weight.allometric.power.sp5 = 3.015816743
species.length2weight.allometric.power.sp6 = 3.256967117
species.length2weight.allometric.power.sp7 = 3.240604905
species.length2weight.allometric.power.sp8 = 3.206546722
species.length2weight.allometric.power.sp9 = 3.242177685
species.length2weight.allometric.power.sp10 = 2.696300975
species.length2weight.allometric.power.sp11 = 3.260833742
species.length2weight.allometric.power.sp12 = 3.20630805
species.length2weight.allometric.power.sp13 = 3.153453379
species.length2weight.allometric.power.sp14 = 3.155472624
species.length2weight.allometric.power.sp15 = 2.806

species.length2weight.condition.factor.sp0 = 0.002000909
species.length2weight.condition.factor.sp1 = 0.006802818
species.length2weight.condition.factor.sp2 = 0.001
species.length2weight.condition.factor.sp3 = 0.0049436
species.length2weight.condition.factor.sp4 = 0.005631986
species.length2weight.condition.factor.sp5 = 0.009380555
species.length2weight.condition.factor.sp6 = 0.00411771
species.length2weight.condition.factor.sp7 = 0.003459801
species.length2weight.condition.factor.sp8 = 0.00451739
species.length2weight.condition.factor.sp9 = 0.004152718
species.length2weight.condition.factor.sp10 = 0.026740884
species.length2weight.condition.factor.sp11 = 0.00349893
species.length2weight.condition.factor.sp12 = 0.005572964
species.length2weight.condition.factor.sp13 = 0.005076828
species.length2weight.condition.factor.sp14 = 0.003626992
species.length2weight.condition.factor.sp15 = 0.0102

species.lInf.sp0 = 30.5
species.lInf.sp1 = 33.6
species.lInf.sp2 = 18.5
species.lInf.sp3 = 15
species.lInf.sp4 = 22.6
species.lInf.sp5 = 48
species.lInf.sp6 = 40
species.lInf.sp7 = 87.7
species.lInf.sp8 = 116
species.lInf.sp9 = 61.4
species.lInf.sp10 = 35
species.lInf.sp11 = 38.1
species.lInf.sp12 = 33
species.lInf.sp13 = 46
species.lInf.sp14 = 127.5
species.lInf.sp15 = 18

species.maturity.size.sp0 = 22.1
species.maturity.size.sp1 = 30
species.maturity.size.sp2 = 13.5
species.maturity.size.sp3 = 11.5
species.maturity.size.sp4 = 15
species.maturity.size.sp5 = 22.1
species.maturity.size.sp6 = 18.8
species.maturity.size.sp7 = 55.4
species.maturity.size.sp8 = 37.8
species.maturity.size.sp9 = 33.5
species.maturity.size.sp10 = 22
species.maturity.size.sp11 = 20.2
species.maturity.size.sp12 = 13
species.maturity.size.sp13 = 25
species.maturity.size.sp14 = 40
species.maturity.size.sp15 = 4.5

species.t0.sp0 = -0.69
species.t0.sp1 = -0.18
species.t0.sp2 = 0
species.t0.sp3 = -0.49
species.t0.sp4 = 0
species.t0.sp5 = -1
species.t0.sp6 = -3
species.t0.sp7 = -0.04
species.t0.sp8 = 0.08
species.t0.sp9 = -0.25
species.t0.sp10 = -1.58
species.t0.sp11 = -0.5
species.t0.sp12 = -0.89
species.t0.sp13 = 0
species.t0.sp14 = -1.13
species.t0.sp15 = -1.14

species.vonbertalanffy.threshold.age.sp0 = 1
species.vonbertalanffy.threshold.age.sp1 = 1
species.vonbertalanffy.threshold.age.sp2 = 1
species.vonbertalanffy.threshold.age.sp3 = 1
species.vonbertalanffy.threshold.age.sp4 = 1
species.vonbertalanffy.threshold.age.sp5 = 1
species.vonbertalanffy.threshold.age.sp6 = 1
species.vonbertalanffy.threshold.age.sp7 = 1
species.vonbertalanffy.threshold.age.sp8 = 1
species.vonbertalanffy.threshold.age.sp9 = 1
species.vonbertalanffy.threshold.age.sp10 = 1
species.vonbertalanffy.threshold.age.sp11 = 1
species.vonbertalanffy.threshold.age.sp12 = 1
species.vonbertalanffy.threshold.age.sp13 = 1
species.vonbertalanffy.threshold.age.sp14 = 1
species.vonbertalanffy.threshold.age.sp15 = 1

species.delta.lmax.factor.sp0 = 1.5
species.delta.lmax.factor.sp1 = 1.5
species.delta.lmax.factor.sp2 = 1.5
species.delta.lmax.factor.sp3 = 1.5
species.delta.lmax.factor.sp4 = 1.5
species.delta.lmax.factor.sp5 = 1.5
species.delta.lmax.factor.sp6 = 1.5
species.delta.lmax.factor.sp7 = 1.5
species.delta.lmax.factor.sp8 = 1.5
species.delta.lmax.factor.sp9 = 1.5
species.delta.lmax.factor.sp10 = 1.5
species.delta.lmax.factor.sp11 = 1.5
species.delta.lmax.factor.sp12 = 1.5
species.delta.lmax.factor.sp13 = 1.5
species.delta.lmax.factor.sp14 = 1.5
species.delta.lmax.factor.sp15 = 1.5

species.type.sp0 = focal
species.type.sp1 = focal
species.type.sp2 = focal
species.type.sp3 = focal
species.type.sp4 = focal
species.type.sp5 = focal
species.type.sp6 = focal
species.type.sp7 = focal
species.type.sp8 = focal
species.type.sp9 = focal
species.type.sp10 = focal
species.type.sp11 = focal
species.type.sp12 = focal
species.type.sp13 = focal
species.type.sp14 = focal
species.type.sp15 = focal

species.type.sp16 = resource
species.type.sp17 = resource
species.type.sp18 = resource
species.type.sp19 = resource
species.type.sp20 = resource
species.type.sp21 = resource
species.type.sp22 = resource
species.type.sp23 = resource
species.type.sp24 = resource
species.type.sp25 = resource

# Reproduction parameters ------------------------------------------------------

reproduction.season.file.sp0 = repro_new_param/reproduction-seasonality-sp0.csv
reproduction.season.file.sp1 = repro_new_param/reproduction-seasonality-sp1.csv
reproduction.season.file.sp2 = repro_new_param/reproduction-seasonality-sp2.csv
reproduction.season.file.sp3 = repro_new_param/reproduction-seasonality-sp3.csv
reproduction.season.file.sp4 = repro_new_param/reproduction-seasonality-sp4.csv
reproduction.season.file.sp5 = repro_new_param/reproduction-seasonality-sp5.csv
reproduction.season.file.sp6 = repro_new_param/reproduction-seasonality-sp6.csv
reproduction.season.file.sp7 = repro_new_param/reproduction-seasonality-sp7.csv
reproduction.season.file.sp8 = repro_new_param/reproduction-seasonality-sp8.csv
reproduction.season.file.sp9 = repro_new_param/reproduction-seasonality-sp9.csv
reproduction.season.file.sp10 = repro_new_param/reproduction-seasonality-sp10.csv
reproduction.season.file.sp11 = repro_new_param/reproduction-seasonality-sp11.csv
reproduction.season.file.sp12 = repro_new_param/reproduction-seasonality-sp12.csv
reproduction.season.file.sp13 = repro_new_param/reproduction-seasonality-sp13.csv
reproduction.season.file.sp14 = repro_new_param/reproduction-seasonality-sp14.csv
reproduction.season.file.sp15 = repro_new_param/reproduction-seasonality-sp15.csv

# Starvation parameters --------------------------------------------------------

mortality.starvation.rate.max.sp0 = 200
mortality.starvation.rate.max.sp1 = 200
mortality.starvation.rate.max.sp2 = 200
mortality.starvation.rate.max.sp3 = 200
mortality.starvation.rate.max.sp4 = 200
mortality.starvation.rate.max.sp5 = 200
mortality.starvation.rate.max.sp6 = 200
mortality.starvation.rate.max.sp7 = 200
mortality.starvation.rate.max.sp8 = 200
mortality.starvation.rate.max.sp9 = 200
mortality.starvation.rate.max.sp10 = 200
mortality.starvation.rate.max.sp11 = 200
mortality.starvation.rate.max.sp12 = 200
mortality.starvation.rate.max.sp13 = 200
mortality.starvation.rate.max.sp14 = 200
mortality.starvation.rate.max.sp15 = 200

# LTL parameters ---------------------------------------------------------------

species.file.sp16 = LTL_ns_corrected_one.nc
species.file.sp17 = LTL_ns_corrected_one.nc
species.file.sp18 = LTL_ns_corrected_one.nc
species.file.sp19 = LTL_ns_corrected_one.nc
species.file.sp20 = LTL_ns_corrected_one.nc
species.file.sp21 = LTL_ns_corrected_one.nc
species.file.sp22 = LTL_ns_corrected_one.nc
species.file.sp23 = LTL_ns_corrected_one.nc
species.file.sp24 = LTL_ns_corrected_one.nc
species.file.sp25 = LTL_ns_corrected_one.nc

#species.accessibility2fish.sp16 = 0
#species.accessibility2fish.sp17 = 0
#species.accessibility2fish.sp18 = 0
#species.accessibility2fish.sp19 = 0
#species.accessibility2fish.sp20 = 0
#species.accessibility2fish.sp21 = 0
#species.accessibility2fish.sp22 = 0
#species.accessibility2fish.sp23 = 0
#species.accessibility2fish.sp24 = 0
#species.accessibility2fish.sp25 = 0

species.name.sp16 = MicroPhytoplankton
species.name.sp17 = Diatoms
species.name.sp18 = HeterotrophicFlagellates
species.name.sp19 = Microzooplankton
species.name.sp20 = Mesozooplankton
species.name.sp21 = SuspensionFeeders
species.name.sp22 = DepositFeeders
species.name.sp23 = Meiobenthos
species.name.sp24 = LargeBenthos
species.name.sp25 = VeryLargeBenthos

species.size.max.sp16 = 0.002
species.size.max.sp17 = 0.02
species.size.max.sp18 = 0.002
species.size.max.sp19 = 0.02
species.size.max.sp20 = 0.5
species.size.max.sp21 = 5
species.size.max.sp22 = 2
species.size.max.sp23 = 0.5
species.size.max.sp24 = 10
species.size.max.sp25 = 15

species.size.min.sp16 = 0.0002
species.size.min.sp17 = 0.001
species.size.min.sp18 = 0.0002
species.size.min.sp19 = 0.002
species.size.min.sp20 = 0.02
species.size.min.sp21 = 2
species.size.min.sp22 = 0.1
species.size.min.sp23 = 0.0045
species.size.min.sp24 = 5
species.size.min.sp25 = 10

species.TL.sp16 = 1
species.TL.sp17 = 1
species.TL.sp18 = 2
species.TL.sp19 = 2
species.TL.sp20 = 2.3
species.TL.sp21 = 2
species.TL.sp22 = 2
species.TL.sp23 = 2
species.TL.sp24 = 2.6
species.TL.sp25 = 3.6

species.biomass.nsteps.year.sp16 = 24
species.biomass.nsteps.year.sp17 = 24
species.biomass.nsteps.year.sp18 = 24
species.biomass.nsteps.year.sp19 = 24
species.biomass.nsteps.year.sp20 = 24
species.biomass.nsteps.year.sp21 = 24
species.biomass.nsteps.year.sp22 = 24
species.biomass.nsteps.year.sp23 = 24
species.biomass.nsteps.year.sp24 = 24
species.biomass.nsteps.year.sp25 = 24

# Predation parameters ---------------------------------------------------------

predation.accessibility.file = predation-accessibility.csv
predation.accessibility.stage.structure = age
predation.accessibility.stage.threshold.sp0 = NULL
predation.accessibility.stage.threshold.sp1 = NULL
predation.accessibility.stage.threshold.sp2 = 0.25
predation.accessibility.stage.threshold.sp3 = NULL
predation.accessibility.stage.threshold.sp4 = 0.25
predation.accessibility.stage.threshold.sp5 = 0.25
predation.accessibility.stage.threshold.sp6 = 0.25
predation.accessibility.stage.threshold.sp7 = 0.25
predation.accessibility.stage.threshold.sp8 = 0.25
predation.accessibility.stage.threshold.sp9 = 0.25
predation.accessibility.stage.threshold.sp10 = NULL
predation.accessibility.stage.threshold.sp11 = 0.25
predation.accessibility.stage.threshold.sp12 = 0.25
predation.accessibility.stage.threshold.sp13 = 0.25
predation.accessibility.stage.threshold.sp14 = 0.25
predation.accessibility.stage.threshold.sp15 = NULL

predation.efficiency.critical.sp0 = 0.57
predation.efficiency.critical.sp1 = 0.57
predation.efficiency.critical.sp10 = 0.57
predation.efficiency.critical.sp11 = 0.57
predation.efficiency.critical.sp12 = 0.57
predation.efficiency.critical.sp13 = 0.57
predation.efficiency.critical.sp2 = 0.57
predation.efficiency.critical.sp3 = 0.57
predation.efficiency.critical.sp4 = 0.57
predation.efficiency.critical.sp5 = 0.57
predation.efficiency.critical.sp6 = 0.57
predation.efficiency.critical.sp7 = 0.57
predation.efficiency.critical.sp8 = 0.57
predation.efficiency.critical.sp9 = 0.57
predation.efficiency.critical.sp14 = 0.57
predation.efficiency.critical.sp15 = 0.57

# old parameter name and value
# predation.ingestion.rate.max.sp0 = 3.5
# predation.ingestion.rate.max.sp1 = 3.5
# predation.ingestion.rate.max.sp10 = 3.5
# predation.ingestion.rate.max.sp11 = 3.5
# predation.ingestion.rate.max.sp12 = 3.5
# predation.ingestion.rate.max.sp13 = 3.5
# predation.ingestion.rate.max.sp2 = 3.5
# predation.ingestion.rate.max.sp3 = 3.5
# predation.ingestion.rate.max.sp4 = 3.5
# predation.ingestion.rate.max.sp5 = 3.5
# predation.ingestion.rate.max.sp6 = 3.5
# predation.ingestion.rate.max.sp7 = 3.5
# predation.ingestion.rate.max.sp8 = 3.5
# predation.ingestion.rate.max.sp9 = 3.5
# predation.ingestion.rate.max.sp14 = 3.5
# predation.ingestion.rate.max.sp15 = 3.5

predation.predPrey.sizeRatio.max.sp0 = 5
predation.predPrey.sizeRatio.max.sp1 = 2.5
predation.predPrey.sizeRatio.max.sp2 = 8
predation.predPrey.sizeRatio.max.sp3 = 8
predation.predPrey.sizeRatio.max.sp4 = 1.5
predation.predPrey.sizeRatio.max.sp5 = 5,5
predation.predPrey.sizeRatio.max.sp6 = 10,10
predation.predPrey.sizeRatio.max.sp7 = 2,2
predation.predPrey.sizeRatio.max.sp8 = 2,1.3
predation.predPrey.sizeRatio.max.sp9 = 1.7,1.7
predation.predPrey.sizeRatio.max.sp10 = 2.5
predation.predPrey.sizeRatio.max.sp11 = 1.5,1.1,1.1
predation.predPrey.sizeRatio.max.sp12 = 5,5
predation.predPrey.sizeRatio.max.sp13 = 2,2
predation.predPrey.sizeRatio.max.sp14 = 2,1.1
predation.predPrey.sizeRatio.max.sp15 = 8

# Test
predation.predPrey.sizeRatio.min.sp0 = 250
predation.predPrey.sizeRatio.min.sp1 = 150
predation.predPrey.sizeRatio.min.sp2 = 1000
predation.predPrey.sizeRatio.min.sp3 = 1000
predation.predPrey.sizeRatio.min.sp4 = 500
predation.predPrey.sizeRatio.min.sp5 = 100,50
predation.predPrey.sizeRatio.min.sp6 = 100,50
predation.predPrey.sizeRatio.min.sp7 = 50,15
predation.predPrey.sizeRatio.min.sp8 = 50,15
predation.predPrey.sizeRatio.min.sp9 = 50,15
predation.predPrey.sizeRatio.min.sp10 = 150
predation.predPrey.sizeRatio.min.sp11 = 100,20,10
predation.predPrey.sizeRatio.min.sp12 = 100,50
predation.predPrey.sizeRatio.min.sp13 = 100,30
predation.predPrey.sizeRatio.min.sp14 = 50,15
predation.predPrey.sizeRatio.min.sp15 = 1000

predation.predPrey.stage.structure = size
predation.predPrey.stage.threshold.sp0 = NULL
predation.predPrey.stage.threshold.sp1 = NULL
predation.predPrey.stage.threshold.sp2 = NULL
predation.predPrey.stage.threshold.sp3 = NULL
predation.predPrey.stage.threshold.sp4 = NULL
predation.predPrey.stage.threshold.sp5 = 8
predation.predPrey.stage.threshold.sp6 = 8
predation.predPrey.stage.threshold.sp7 = 20
predation.predPrey.stage.threshold.sp8 = 12
predation.predPrey.stage.threshold.sp9 = 12
predation.predPrey.stage.threshold.sp10 = NULL
predation.predPrey.stage.threshold.sp11 = 10,40
predation.predPrey.stage.threshold.sp12 = 30
predation.predPrey.stage.threshold.sp13 = 12
predation.predPrey.stage.threshold.sp14 = 12
predation.predPrey.stage.threshold.sp15 = NULL

# Out mortality parameters -----------------------------------------------------

mortality.out.rate.sp0 = 0
mortality.out.rate.sp1 = 0
mortality.out.rate.sp2 = 0
mortality.out.rate.sp3 = 0
mortality.out.rate.sp4 = 0
mortality.out.rate.sp5 = 0
mortality.out.rate.sp6 = 0
mortality.out.rate.sp7 = 0
mortality.out.rate.sp8 = 0
mortality.out.rate.sp9 = 0
mortality.out.rate.sp10 = 0
mortality.out.rate.sp11 = 0
mortality.out.rate.sp12 = 0
mortality.out.rate.sp13 = 0
mortality.out.rate.sp14 = 0

# Movements parameters ---------------------------------------------------------

movement.netcdf.enabled = FALSE
movement.checks.enabled = FALSE

movement.distribution.method.sp0 = maps
movement.distribution.method.sp1 = maps
movement.distribution.method.sp2 = maps
movement.distribution.method.sp3 = maps
movement.distribution.method.sp4 = maps
movement.distribution.method.sp5 = maps
movement.distribution.method.sp6 = maps
movement.distribution.method.sp7 = maps
movement.distribution.method.sp8 = maps
movement.distribution.method.sp9 = maps
movement.distribution.method.sp10 = maps
movement.distribution.method.sp11 = maps
movement.distribution.method.sp12 = maps
movement.distribution.method.sp13 = maps
movement.distribution.method.sp14 = maps
movement.distribution.method.sp15 = maps

movement.lastAge.map0 = 17
movement.initialAge.map0 = 2.9
movement.file.map0 = maps/0Herring_adu.csv
movement.steps.map0 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map0 = AtlanticHerring

movement.lastAge.map1 = 2.9
movement.initialAge.map1 = 0.25
movement.file.map1 = maps/0Herring_juv.csv
movement.steps.map1 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map1 = AtlanticHerring

movement.lastAge.map2 = 17
movement.initialAge.map2 = 2.45
movement.file.map2 = maps/1Mackerel_adu.csv
movement.steps.map2 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map2 = AtlanticMackerel

movement.lastAge.map3 = 2.45
movement.initialAge.map3 = 0.25
movement.file.map3 = maps/1Mackerel_juv.csv
movement.steps.map3 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map3 = AtlanticMackerel

movement.lastAge.map4 = 10
movement.initialAge.map4 = 2.16
movement.file.map4 = maps/2Sandeel_adu.csv
movement.steps.map4 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map4 = Sandeel

movement.lastAge.map5 = 2.16
movement.initialAge.map5 = 0
movement.file.map5 = maps/2Sandeel_juv.csv
movement.steps.map5 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map5 = Sandeel

movement.lastAge.map6 = 10
movement.initialAge.map6 = 2.28
movement.file.map6 = maps/3EuropeanSprat_adu.csv
movement.steps.map6 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map6 = EuropeanSprat

movement.lastAge.map7 = 2.28
movement.initialAge.map7 = 0.25
movement.file.map7 = maps/3EuropeanSprat_juv.csv
movement.steps.map7 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map7 = EuropeanSprat

movement.lastAge.map8 = 7
movement.initialAge.map8 = 1.97
movement.file.map8 = maps/4NorwayPout_adu.csv
movement.steps.map8 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map8 = NorwayPout

movement.lastAge.map9 = 1.97
movement.initialAge.map9 = 0.25
movement.file.map9 = maps/4NorwayPout_juv.csv
movement.steps.map9 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map9 = NorwayPout

movement.lastAge.map10 = 25
movement.initialAge.map10 = 3.82
movement.file.map10 = maps/5EuropeanPlaice_adu.csv
movement.steps.map10 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map10 = EuropeanPlaice

movement.lastAge.map11 = 3.82
movement.initialAge.map11 = 0.25
movement.file.map11 = maps/5EuropeanPlaice_juv.csv
movement.steps.map11 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map11 = EuropeanPlaice

movement.lastAge.map12 = 26
movement.initialAge.map12 = 2.98
movement.file.map12 = maps/6CommonSole_adu.csv
movement.steps.map12 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map12 = CommonSole

movement.lastAge.map13 = 2.98
movement.initialAge.map13 = 0.25
movement.file.map13 = maps/6CommonSole_juv.csv
movement.steps.map13 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map13 = CommonSole

movement.lastAge.map14 = 25
movement.initialAge.map14 = 5.79
movement.file.map14 = maps/7Saithe_adu.csv
movement.steps.map14 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map14 = Saithe

movement.lastAge.map15 = 5.79
movement.initialAge.map15 = 0.25
movement.file.map15 = maps/7Saithe_juv.csv
movement.steps.map15 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map15 = Saithe

movement.lastAge.map16 = 25
movement.initialAge.map16 = 3.9
movement.file.map16 = maps/8AtlanticCod_adu.csv
movement.steps.map16 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map16 = AtlanticCod

movement.lastAge.map17 = 3.9
movement.initialAge.map17 = 0.25
movement.file.map17 = maps/8AtlanticCod_juv.csv
movement.steps.map17 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map17 = AtlanticCod

movement.lastAge.map18 = 20
movement.initialAge.map18 = 2.6
movement.file.map18 = maps/9Haddock_adu.csv
movement.steps.map18 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map18 = Haddock

movement.lastAge.map19 = 2.6
movement.initialAge.map19 = 0.25
movement.file.map19 = maps/9Haddock_juv.csv
movement.steps.map19 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map19 = Haddock

movement.lastAge.map20 = 15
movement.initialAge.map20 = 1.76
movement.file.map20 = maps/10HorseMackerel_adu.csv
movement.steps.map20 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map20 = HorseMackerel

movement.lastAge.map21 = 1.76
movement.initialAge.map21 = 0
movement.file.map21 = maps/10HorseMackerel_juv.csv
movement.steps.map21 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map21 = HorseMackerel

movement.lastAge.map22 = 20
movement.initialAge.map22 = 0.25
movement.file.map22 = maps/11Whiting_adu.csv
movement.steps.map22 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map22 = Whiting

movement.lastAge.map23 = 12
movement.initialAge.map23 = 0.25
movement.file.map23 = maps/12CommonDab_adu.csv
movement.steps.map23 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map23 = CommonDab

movement.lastAge.map24 = 21
movement.initialAge.map24 = 0
movement.file.map24 = maps/13GreyGurnard_adu.csv
movement.steps.map24 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map24 = GreyGurnard

movement.lastAge.map25 = 20
movement.initialAge.map25 = 3.4
movement.file.map25 = maps/14Hake_adu.csv
movement.steps.map25 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map25 = Hake

movement.lastAge.map26 = 3.4
movement.initialAge.map26 = 0
movement.file.map26 = maps/14Hake_juv.csv
movement.steps.map26 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map26 = Hake

movement.lastAge.map27 = 3
movement.initialAge.map27 = 0
movement.file.map27 = maps/15Shrimp_adu.csv
movement.steps.map27 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map27 = Shrimp

movement.lastAge.map28 = 0.25
movement.initialAge.map28 = 0
movement.file.map28 = maps/0Herring_egg.csv
movement.steps.map28 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map28 = AtlanticHerring

movement.lastAge.map29 = 0.25
movement.initialAge.map29 = 0
movement.file.map29 = maps/1Mackerel_egg.csv
movement.steps.map29 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map29 = AtlanticMackerel

movement.lastAge.map30 = 0.25
movement.initialAge.map30 = 0
movement.file.map30 = maps/3EuropeanSprat_egg.csv
movement.steps.map30 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map30 = EuropeanSprat

movement.lastAge.map31 = 0.25
movement.initialAge.map31 = 0
movement.file.map31 = maps/4NorwayPout_egg.csv
movement.steps.map31 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map31 = NorwayPout

movement.lastAge.map32 = 0.25
movement.initialAge.map32 = 0
movement.file.map32 = maps/5EuropeanPlaice_egg.csv
movement.steps.map32 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map32 = EuropeanPlaice

movement.lastAge.map33 = 0.25
movement.initialAge.map33 = 0
movement.file.map33 = maps/6CommonSole_egg.csv
movement.steps.map33 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map33 = CommonSole

movement.lastAge.map34 = 0.25
movement.initialAge.map34 = 0
movement.file.map34 = maps/7Saithe_egg.csv
movement.steps.map34 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map34 = Saithe

movement.lastAge.map35 = 0.25
movement.initialAge.map35 = 0
movement.file.map35 = maps/8AtlanticCod_egg.csv
movement.steps.map35 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map35 = AtlanticCod

movement.lastAge.map36 = 0.25
movement.initialAge.map36 = 0
movement.file.map36 = maps/9Haddock_egg.csv
movement.steps.map36 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map36 = Haddock

movement.lastAge.map37 = 0.25
movement.initialAge.map37 = 0
movement.file.map37 = maps/10HorseMackerel_egg.csv
movement.steps.map37 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map37 = HorseMackerel

movement.lastAge.map38 = 0.25
movement.initialAge.map38 = 0
movement.file.map38 = maps/11Whiting_egg.csv
movement.steps.map38 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map38 = Whiting

movement.lastAge.map39 = 0.25
movement.initialAge.map39 = 0
movement.file.map39 = maps/12CommonDab_egg.csv
movement.steps.map39 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map39 = CommonDab

movement.lastAge.map40 = 0.25
movement.initialAge.map40 = 0
movement.file.map40 = maps/14Hake_egg.csv
movement.steps.map40 = 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
movement.species.map40 = Hake

movement.randomwalk.range.sp6 = 1
movement.randomwalk.range.sp7 = 1
movement.randomwalk.range.sp8 = 1
movement.randomwalk.range.sp9 = 1
movement.randomwalk.range.sp10 = 1
movement.randomwalk.range.sp11 = 1
movement.randomwalk.range.sp12 = 1
movement.randomwalk.range.sp13 = 1
movement.randomwalk.range.sp14 = 1
movement.randomwalk.range.sp15 = 1

# OLD Fishing parameters -------------------------------------------------------

mortality.fishing.recruitment.size.sp0 = 16.5
mortality.fishing.recruitment.size.sp1 = 21
mortality.fishing.recruitment.size.sp2 = 12
mortality.fishing.recruitment.size.sp3 = 10
mortality.fishing.recruitment.size.sp4 = 12
mortality.fishing.recruitment.size.sp5 = 15
mortality.fishing.recruitment.size.sp6 = 12
mortality.fishing.recruitment.size.sp7 = 29.5
mortality.fishing.recruitment.size.sp8 = 21
mortality.fishing.recruitment.size.sp9 = 18
mortality.fishing.recruitment.size.sp10 = 18
mortality.fishing.recruitment.size.sp11 = 16
mortality.fishing.recruitment.size.sp12 = 12
mortality.fishing.recruitment.size.sp13 = 15
mortality.fishing.recruitment.size.sp14 = 23
mortality.fishing.recruitment.size.sp15 = 5

mortality.fishing.rate.byDt.bySize.file.sp0 = fishing_rate/mortality.fishing.rate.byDt.bySize.file.sp0.csv
mortality.fishing.rate.byDt.bySize.file.sp1 = F_Fstart1_Fmax0.5_sp1.csv
mortality.fishing.rate.byDt.bySize.file.sp2 = fishing_rate/mortality.fishing.rate.byDt.bySize.file.sp2.csv
mortality.fishing.rate.byDt.bySize.file.sp3 = fishing_rate/mortality.fishing.rate.byDt.bySize.file.sp3.csv
mortality.fishing.rate.byDt.bySize.file.sp4 = fishing_rate/mortality.fishing.rate.byDt.bySize.file.sp4.csv
mortality.fishing.rate.byDt.bySize.file.sp5 = fishing_rate/mortality.fishing.rate.byDt.bySize.file.sp5.csv
mortality.fishing.rate.byDt.bySize.file.sp6 = fishing_rate/mortality.fishing.rate.byDt.bySize.file.sp6.csv
mortality.fishing.rate.byDt.bySize.file.sp7 = fishing_rate/mortality.fishing.rate.byDt.bySize.file.sp7.csv
mortality.fishing.rate.byDt.bySize.file.sp8 = fishing_rate/mortality.fishing.rate.byDt.bySize.file.sp8.csv
mortality.fishing.rate.byDt.bySize.file.sp9 = fishing_rate/mortality.fishing.rate.byDt.bySize.file.sp9.csv
mortality.fishing.rate.byDt.bySize.file.sp10 = fishing_rate/mortality.fishing.rate.byDt.bySize.file.sp10.csv
mortality.fishing.rate.byDt.bySize.file.sp11 = fishing_rate/mortality.fishing.rate.byDt.bySize.file.sp11.csv
mortality.fishing.rate.byDt.bySize.file.sp12 = fishing_rate/mortality.fishing.rate.byDt.bySize.file.sp12.csv
mortality.fishing.rate.byDt.bySize.file.sp13 = fishing_rate/mortality.fishing.rate.byDt.bySize.file.sp13.csv
mortality.fishing.rate.byDt.bySize.file.sp14 = fishing_rate/mortality.fishing.rate.byDt.bySize.file.sp14.csv
mortality.fishing.rate.byDt.bySize.file.sp15 = fishing_rate/mortality.fishing.rate.byDt.bySize.file.sp15.csv

# Init pop parameters ----------------------------------------------------------

# Osmose 3 Update 2 (2015/03/01) Deprecated parameter population.initialization.method
# population.initialization.method,biomass
population.seeding.biomass.sp0 = 4251786.889
population.seeding.biomass.sp1 = 1719253.293
population.seeding.biomass.sp2 = 1734689.272
population.seeding.biomass.sp3 = 1309145.875
population.seeding.biomass.sp4 = 393338.7041
population.seeding.biomass.sp5 = 969445.942
population.seeding.biomass.sp6 = 78223.22565
population.seeding.biomass.sp7 = 354496.1885
population.seeding.biomass.sp8 = 205978.603
population.seeding.biomass.sp9 = 379206.8846
population.seeding.biomass.sp10 = 330000
population.seeding.biomass.sp11 = 362789.125
population.seeding.biomass.sp12 = 201125.2721
population.seeding.biomass.sp13 = 43900
population.seeding.biomass.sp14 = 61687.2426
population.seeding.biomass.sp15 = 125904



# Simulation parameters --------------------------------------------------------

simulation.ncpu = 1
simulation.nsimulation = 20
simulation.nresource = 10
simulation.nspecies = 16
simulation.time.ndtperyear = 24
# simulation.time.nstep = 24
#simulation.time.nyear,1
# 160 y for python
simulation.time.nyear = 200
#simulation.nschool,30
simulation.nschool.sp0 = 25
simulation.nschool.sp1 = 38
simulation.nschool.sp2 = 38
simulation.nschool.sp3 = 38
simulation.nschool.sp4 = 30
simulation.nschool.sp5 = 50
simulation.nschool.sp6 = 50
simulation.nschool.sp7 = 30
simulation.nschool.sp8 = 25
simulation.nschool.sp9 = 38
simulation.nschool.sp10 = 38
simulation.nschool.sp11 = 38
simulation.nschool.sp12 = 50
simulation.nschool.sp13 = 30
simulation.nschool.sp14 = 20
simulation.nschool.sp15 = 25



# Model initialisation ----------------------------------------------------

population.seeding.year.max = 50
population.initialization.relativebiomass.enabled = FALSE
# 
# # this file provide the parameters to generate the initial conditions
# osmose.configuration.initialisation.setup = input/initialisation_setup_parameters.osm
# # this file provide the actual initial conditions
# osmose.configuration.initial.conditions = input/initial_conditions.osm


# Output parameters ------------------------------------------------------------

# Main output parameters
output.start.year = 0
output.file.prefix = ns
# 1 means records every time step (24), 24 means it records once a year
output.recordfrequency.ndt = 24
output.netcdf.format=netcdf3
# CSV separator (COMA, SEMICOLON, EQUALS, COLON, TAB
#output.csv.separator = COMA

# Save restart file
# Osmose 4.4.0 - Renamed parameter output.restart.enabled into simulation.restart.enabled
output.restart.enabled = TRUE
simulation.restart.enabled = TRUE
# Osmose 4.4.0 - Renamed parameter output.restart.spinup into simulation.restart.spinup.nyear
simulation.restart.spinup.nyear = 199
# Osmose 4.4.0 - Renamed parameter output.restart.recordFrequency.ndt into simulation.restart.recordFrequency.ndt
simulation.restart.recordFrequency.ndt = 24
simulation.restart.file = RESTART/ns_snapshot_step1919.nc.3

# Biomass
output.biomass.enabled = TRUE
output.biomass.bysize.enabled = TRUE
output.biomass.byage.enabled = TRUE
# output.abundance.byweight.enabled = TRUE
# output.biomass.byweight.enabled = TRUE
output.ssb.enabled = TRUE

# Abundance
output.abundance.enabled = TRUE
output.abundance.bysize.enabled = TRUE
output.abundance.byage.enabled = TRUE
# output.abundance.age1.enabled = FALSE
output.mature.abundance.byage.enabled = TRUE
output.mature.biomass.byage.enabled = TRUE

# Mortality
output.mortality.enabled = TRUE
output.mortality.perSpecies.byAge.enabled = TRUE
output.mortality.perSpecies.bySize.enabled = TRUE
output.mortality.perSpecies.byage.enabled = TRUE
output.mortality.perSpecies.bysize.enabled = TRUE
output.mortality.additional.bySize.enabled = TRUE
output.mortality.additional.byAge.enabled = TRUE
# output.mortality.additionalN.bySize.enabled = FALSE
# output.mortality.additionalN.byAge.enabled = FALSE

# Yield
output.yield.biomass.enabled = TRUE
output.yield.abundance.enabled = TRUE
# output.yieldN.bySize.enabled = FALSE
output.yield.bySize.enabled = TRUE
# output.yieldN.byAge.enabled = FALSE
output.yield.byAge.enabled = TRUE
# output.yield.abundance.netcdf.enabled = TRUE
# output.yield.biomass.bySize.netcdf.enabled = FALSE
# output.yield.biomass.byage.netcdf.enabled = FALSE
# output.yield.abundance.bySize.netcdf.enabled = FALSE
# output.yield.abundance.byage.netcdf.enabled = FALSE
# output.spatial.yield.biomass.enabled = FALSE
# output.spatial.yield.abundance.enabled = FALSE
output.yield.abundance.bySize.enabled = TRUE
output.yield.biomass.bySize.enabled = TRUE
output.yield.abundance.byAge.enabled = TRUE
output.yield.biomass.byAge.enabled = TRUE

# Size
# output.size.enabled = FALSE
# output.size.catch.enabled = FALSE
output.meanSize.byAge.enabled = TRUE
output.meanWeight.byAge.enabled = TRUE

# TL
# output.TL.enabled = FALSE
# output.TL.catch.enabled = FALSE
# output.biomass.byTL.enabled = FALSE
# output.meanTL.bySize.enabled = FALSE
#output.meanTL.byAge.enabled = FALSE
output.abundance.bytl.enabled = TRUE

# Predation
# output.diet.composition.enabled = FALSE
# output.diet.composition.byAge.enabled = FALSE
# output.diet.composition.bySize.enabled = FALSE
# output.diet.pressure.enabled = FALSE
# output.diet.pressure.byAge.enabled = FALSE
# output.diet.pressure.bySize.enabled = FALSE

# Spatial
# output.spatial.enabled = FALSE
# output.spatial.ltl.enabled = FALSE

# Advanced parameters
output.nschool.enabled = TRUE
output.nschool.byage.enabled = TRUE
# Whether to include step 0 of the simulation in the outputs
output.step0.include = FALSE
# Cutoff for biomass, abundance, mean size and mean trophic level outputs
output.cutoff.enabled = FALSE
output.cutoff.age.sp0 = 1
output.cutoff.age.sp1 = 1
output.cutoff.age.sp2 = 1
output.cutoff.age.sp3 = 1
output.cutoff.age.sp4 = 1
output.cutoff.age.sp5 = 1
output.cutoff.age.sp6 = 1
output.cutoff.age.sp7 = 1
output.cutoff.age.sp8 = 1
output.cutoff.age.sp9 = 1
output.cutoff.age.sp10 = 1
output.cutoff.age.sp11 = 1
output.cutoff.age.sp12 = 1
output.cutoff.age.sp13 = 1
output.cutoff.age.sp14 = 1
output.cutoff.age.sp15 = 1

# Diet output stages
output.diet.success.enabled = TRUE
output.diet.stage.structure = size
output.diet.stage.threshold.sp0 = 10,20,30,40,50
output.diet.stage.threshold.sp1 = 10,20,30,40,50
output.diet.stage.threshold.sp2 = 5,10,20
output.diet.stage.threshold.sp3 = 5,10,20
output.diet.stage.threshold.sp4 = 10,20,30
output.diet.stage.threshold.sp5 = 10,20,30,40,50
output.diet.stage.threshold.sp6 = 10,20,30,40,50
output.diet.stage.threshold.sp7 = 10,20,30,40,50,60,70,80,90,100
output.diet.stage.threshold.sp8 = 10,20,30,40,50,60,70,80,90,100
output.diet.stage.threshold.sp9 = 10,20,30,40,50,60,70,80,90,100
output.diet.stage.threshold.sp10 = 10,20,30,40
output.diet.stage.threshold.sp11 = 10,20,30,40,50,60
output.diet.stage.threshold.sp12 = 10,20,30,40
output.diet.stage.threshold.sp13 = 10,20,30,40,50
output.diet.stage.threshold.sp14 = 10,20,30,40,50,60,70,80,90,100
output.diet.stage.threshold.sp15 = 3,6

# Age tion (year)
output.distrib.byAge.max = 30
output.distrib.byAge.incr = 1
output.distrib.byAge.min = 0

# Size tion (centimetre)
output.distrib.bySize.max = 205
output.distrib.bySize.incr = 10

#netcdf
output.use.netcdf = TRUE
# output.abundance.netcdf.enabled = FALSE
# output.biomass.netcdf.enabled = FALSE
# output.yield.biomass.netcdf.enabled = FALSE
# output.yieldN.netcdf.enabled = FALSE
# output.biomass.bysize.netcdf.enabled = FALSE
# output.biomass.bytl.netcdf.enabled = FALSE
# output.biomass.byage.netcdf.enabled = FALSE
# output.abundance.bysize.netcdf.enabled = FALSE
# output.abundance.byage.netcdf.enabled = FALSE
# output.diet.pressure.netcdf.enabled = FALSE
# output.diet.composition.netcdf.enabled = FALSE
# output.mortality.perSpecies.byage.netcdf.enabled = FALSE
# output.diet.composition.byage.netcdf.enabled = FALSE
# output.diet.composition.bysize.netcdf.enabled = FALSE
# output.yield.bySize.netcdf.enabled = FALSE
# output.yield.byage.netcdf.enabled = FALSE
# output.yieldN.bySize.netcdf.enabled = FALSE
# output.yieldN.byage.netcdf.enabled = FALSE
# output.meanSize.byAge.netcdf.enabled = FALSE

# output.mortality.natural.bySize.enabled = FALSE
# output.mortality.natural.byAge.enabled = FALSE
# output.mortality.naturalN.bySize.enabled = FALSE
# output.mortality.naturalN.byAge.enabled = FALSE
# output.yieldN.enabled = FALSE
# output.weight.enabled = FALSE

#Spatial
# output.spatial.biomass.enabled = FALSE
# output.spatial.abundance.enabled = TRUE
# output.spatial.yield.enabled = FALSE
# output.spatial.yieldN.enabled = FALSE
# output.spatial.size.enabled = FALSE
# output.spatialtl.enabled = FALSE
# output.spatialsizespecies.enabled = FALSE
# output.spatialagespecies.enabled = FALSE

#Individual
# output.individual.enabled = FALSE

# # Fisheries
# output.fishing.accessible.biomass = FALSE # was TRUE
# output.fishing.harvested.biomass  = FALSE # was TRUE
# output.fisheries.enabled = FALSE # was TRUE
# output.fisheries.byweight.enabled = FALSE  # was TRUE
# output.fisheries.bysize.enabled = FALSE
# output.fisheries.byage.enabled = FALSE
# 
# # Economic module
# output.economic.harvesting.costs.enabled = FALSE # was TRUE
# output.economic.fish.prices.enabled = FALSE # was TRUE
# output.economic.fisherman.profit.enabled = FALSE # was TRUE
# output.economic.profit.margin.enabled = FALSE # was TRUE
# output.economic.utility.enabled = FALSE # was TRUE

# Bioen parameters -------------------------------------------------------------

# Osmose 4.4.0 - Renamed parameter simulation.bioen.enabled into module.bioenergetics.enabled
module.bioenergetics.enabled = TRUE

species.bioenergetics.model.sp0 = full
species.bioenergetics.model.sp1 = full
species.bioenergetics.model.sp2 = full
species.bioenergetics.model.sp3 = full
species.bioenergetics.model.sp4 = full
species.bioenergetics.model.sp5 = full
species.bioenergetics.model.sp6 = full
species.bioenergetics.model.sp7 = full
species.bioenergetics.model.sp8 = full
species.bioenergetics.model.sp9 = full
species.bioenergetics.model.sp10 = full
species.bioenergetics.model.sp11 = full
species.bioenergetics.model.sp12 = full
species.bioenergetics.model.sp13 = full
species.bioenergetics.model.sp14 = full
species.bioenergetics.model.sp15 = full

output.fecundity.byage.enabled = TRUE
output.age.death = TRUE
# output.bioen.mature.size.enabled = FALSE
output.bioen.mature.age.enabled = TRUE
# output.bioen.meanGonadweight.enabled = FALSE
# output.bioen.sizeInf.enabled = FALSE
# output.bioen.ingest.enabled = FALSE
# output.bioen.maint.enabled = FALSE
# output.bioen.enet.enabled = FALSE
# output.bioen.rho.enabled = FALSE

output.ingest.byAge.enabled = TRUE
# output.ingest.bySize.enabled = FALSE
output.maintenance.byAge.enabled = TRUE
# output.maintenance.bySize.enabled = FALSE
output.enet.byAge.enabled = TRUE
# output.enet.bySize.enabled = FALSE
output.meanGonadWeight.byAge.enabled = TRUE
output.meanSomaticWeight.byAge.enabled = TRUE
output.rho.byAge.enabled = TRUE
# output.rho.bySize.enabled = FALSE


###### oxygen parameters for the fO2 function
species.oxygen.c1.sp0 = 1
species.oxygen.c1.sp1 = 1
species.oxygen.c1.sp2 = 1
species.oxygen.c1.sp3 = 1
species.oxygen.c1.sp4 = 1.48
species.oxygen.c1.sp5 = 2
species.oxygen.c1.sp6 = 2.6
species.oxygen.c1.sp7 = 1.48
species.oxygen.c1.sp8 = 1.48
species.oxygen.c1.sp9 = 1.48
species.oxygen.c1.sp10 = 1
species.oxygen.c1.sp11 = 1.48
species.oxygen.c1.sp12 = 2.3095
species.oxygen.c1.sp13 = 1.48
species.oxygen.c1.sp14 = 1.48
species.oxygen.c1.sp15 = 1.475

species.oxygen.c2.sp0 = 8.1
species.oxygen.c2.sp1 = 8.1
species.oxygen.c2.sp2 = 7.3
species.oxygen.c2.sp3 = 8.1
species.oxygen.c2.sp4 = 54.6
species.oxygen.c2.sp5 = 99
species.oxygen.c2.sp6 = 179.4
species.oxygen.c2.sp7 = 54.6
species.oxygen.c2.sp8 = 54.6
species.oxygen.c2.sp9 = 54.6
species.oxygen.c2.sp10 = 8.1
species.oxygen.c2.sp11 = 54.6
species.oxygen.c2.sp12 = 139.189
species.oxygen.c2.sp13 = 54.6
species.oxygen.c2.sp14 = 54.6
species.oxygen.c2.sp15 = 47

############### sets the index of the data file to read.
# 0 = surface," 1 = mean surface-bottom"
species.zlayer.sp0 = 0
species.zlayer.sp1 = 0
species.zlayer.sp2 = 1
species.zlayer.sp3 = 0
species.zlayer.sp4 = 1
species.zlayer.sp5 = 2
species.zlayer.sp6 = 2
species.zlayer.sp7 = 1
species.zlayer.sp8 = 1
species.zlayer.sp9 = 1
species.zlayer.sp10 = 0
species.zlayer.sp11 = 1
species.zlayer.sp12 = 2
species.zlayer.sp13 = 1
species.zlayer.sp14 = 1
species.zlayer.sp15 = 1

### temperature coefficients
species.bioen.mobilized.e.mobi.sp0 = 1.929569294
species.bioen.mobilized.e.mobi.sp1 = 1.674103935
species.bioen.mobilized.e.mobi.sp2 = 1.955097947
species.bioen.mobilized.e.mobi.sp3 = 1.151577457
species.bioen.mobilized.e.mobi.sp4 = 3.20728334
species.bioen.mobilized.e.mobi.sp5 = 1.796432348
species.bioen.mobilized.e.mobi.sp6 = 2.544341599
species.bioen.mobilized.e.mobi.sp7 = 2.640416684
species.bioen.mobilized.e.mobi.sp8 = 2.265304514
species.bioen.mobilized.e.mobi.sp9 = 2.434481128
species.bioen.mobilized.e.mobi.sp10 = 2.279256748
species.bioen.mobilized.e.mobi.sp11 = 1.925563754
species.bioen.mobilized.e.mobi.sp12 = 2.550248837
species.bioen.mobilized.e.mobi.sp13 = 1.979443567
species.bioen.mobilized.e.mobi.sp14 = 2.922773121
species.bioen.mobilized.e.mobi.sp15 = 2.250380933

species.bioen.mobilized.e.D.sp0 = 1.929670119
species.bioen.mobilized.e.D.sp1 = 1.674205724
species.bioen.mobilized.e.D.sp2 = 1.955198854
species.bioen.mobilized.e.D.sp3 = 1.151683864
species.bioen.mobilized.e.D.sp4 = 3.207377613
species.bioen.mobilized.e.D.sp5 = 1.79653389
species.bioen.mobilized.e.D.sp6 = 2.544440319
species.bioen.mobilized.e.D.sp7 = 2.64051254
species.bioen.mobilized.e.D.sp8 = 2.265402914
species.bioen.mobilized.e.D.sp9 = 2.434578036
species.bioen.mobilized.e.D.sp10 = 2.27935564
species.bioen.mobilized.e.D.sp11 = 1.9256644
species.bioen.mobilized.e.D.sp12 = 2.550347483
species.bioen.mobilized.e.D.sp13 = 1.97954404
species.bioen.mobilized.e.D.sp14 = 2.922869154
species.bioen.mobilized.e.D.sp15 = 2.25048059

# Tp in Celcius
species.bioen.mobilized.Tp.sp0 = 43.09747476
species.bioen.mobilized.Tp.sp1 = 46.84507211
species.bioen.mobilized.Tp.sp2 = 43.29032659
species.bioen.mobilized.Tp.sp3 = 68.31343394
species.bioen.mobilized.Tp.sp4 = 23.79511454
species.bioen.mobilized.Tp.sp5 = 46.58860357
species.bioen.mobilized.Tp.sp6 = 35.12117599
species.bioen.mobilized.Tp.sp7 = 27.86946487
species.bioen.mobilized.Tp.sp8 = 34.77337736
species.bioen.mobilized.Tp.sp9 = 30.67261219
species.bioen.mobilized.Tp.sp10 = 39.38575479
species.bioen.mobilized.Tp.sp11 = 41.90914441
species.bioen.mobilized.Tp.sp12 = 34.91749722
species.bioen.mobilized.Tp.sp13 = 41.48788084
species.bioen.mobilized.Tp.sp14 = 30.05714505
species.bioen.mobilized.Tp.sp15 = 38.53796895

species.bioen.maint.e.maint.sp0 = 0.595682223
species.bioen.maint.e.maint.sp1 = 0.333531516
species.bioen.maint.e.maint.sp2 = 0.47339626
species.bioen.maint.e.maint.sp3 = 0.268933469
species.bioen.maint.e.maint.sp4 = 0.173632518
species.bioen.maint.e.maint.sp5 = 0.37530696
species.bioen.maint.e.maint.sp6 = 0.417973605
species.bioen.maint.e.maint.sp7 = 0.2588154
species.bioen.maint.e.maint.sp8 = 0.477614462
species.bioen.maint.e.maint.sp9 = 0.307223157
species.bioen.maint.e.maint.sp10 = 0.246413775
species.bioen.maint.e.maint.sp11 = 0.298152628
species.bioen.maint.e.maint.sp12 = 0.458920549
species.bioen.maint.e.maint.sp13 = 0.308743327
species.bioen.maint.e.maint.sp14 = 0.248132689
species.bioen.maint.e.maint.sp15 = 0.404345929

species.bioen.maint.energy.c_m.sp0 = 1.49e+11
species.bioen.maint.energy.c_m.sp1 = 5034832.882
species.bioen.maint.energy.c_m.sp2 = 540388517.7
species.bioen.maint.energy.c_m.sp3 = 256579.0172
species.bioen.maint.energy.c_m.sp4 = 5816.281889
species.bioen.maint.energy.c_m.sp5 = 15368755.3
species.bioen.maint.energy.c_m.sp6 = 52809960.9
species.bioen.maint.energy.c_m.sp7 = 274469.8036
species.bioen.maint.energy.c_m.sp8 = 2066760649
species.bioen.maint.energy.c_m.sp9 = 1961447.834
species.bioen.maint.energy.c_m.sp10 = 103677.1587
species.bioen.maint.energy.c_m.sp11 = 1357998.53
species.bioen.maint.energy.c_m.sp12 = 286079471.6
species.bioen.maint.energy.c_m.sp13 = 1474805.094
species.bioen.maint.energy.c_m.sp14 = 174241.9494
species.bioen.maint.energy.c_m.sp15 = 31010952.85

################# oxydation mortality
# ???? not yet filled
species.bioen.forage.k_for.sp0 = 0.05
species.bioen.forage.k_for.sp1 = 0.05
species.bioen.forage.k_for.sp2 = 0.05
species.bioen.forage.k_for.sp3 = 0.05
species.bioen.forage.k_for.sp4 = 0.05
species.bioen.forage.k_for.sp5 = 0.05
species.bioen.forage.k_for.sp6 = 0.05
species.bioen.forage.k_for.sp7 = 0.05
species.bioen.forage.k_for.sp8 = 0.05
species.bioen.forage.k_for.sp9 = 0.05
species.bioen.forage.k_for.sp10 = 0.05
species.bioen.forage.k_for.sp11 = 0.05
species.bioen.forage.k_for.sp12 = 0.05
species.bioen.forage.k_for.sp13 = 0.05
species.bioen.forage.k_for.sp14 = 0.05
species.bioen.forage.k_for.sp15 = 0.05

species.bioen.forage.k1_for.sp0 = 0.05395111
species.bioen.forage.k1_for.sp1 = 0.055338
species.bioen.forage.k1_for.sp2 = 0.05219479
species.bioen.forage.k1_for.sp3 = 0.05693232
species.bioen.forage.k1_for.sp4 = 0.05413156
species.bioen.forage.k1_for.sp5 = 0.04972486
species.bioen.forage.k1_for.sp6 = 0.03994102
species.bioen.forage.k1_for.sp7 = 0.05448794
species.bioen.forage.k1_for.sp8 = 0.03114422
species.bioen.forage.k1_for.sp9 = 0.04046042
species.bioen.forage.k1_for.sp10 = 0.03606126
species.bioen.forage.k1_for.sp11 = 0.04107436
species.bioen.forage.k1_for.sp12 = 0.05196456
species.bioen.forage.k1_for.sp13 = 0.04674665
species.bioen.forage.k1_for.sp14 = 0.03922959
species.bioen.forage.k1_for.sp15 = 0.05

species.bioen.forage.k2_for.sp0 = 0.95
species.bioen.forage.k2_for.sp1 = 0.7
species.bioen.forage.k2_for.sp2 = 0.95
species.bioen.forage.k2_for.sp3 = 1.1
species.bioen.forage.k2_for.sp4 = 1.1
species.bioen.forage.k2_for.sp5 = 1.1
species.bioen.forage.k2_for.sp6 = 1.1
species.bioen.forage.k2_for.sp7 = 0.95
species.bioen.forage.k2_for.sp8 = 0.95
species.bioen.forage.k2_for.sp9 = 0.95
species.bioen.forage.k2_for.sp10 = 0.95
species.bioen.forage.k2_for.sp11 = 0.95
species.bioen.forage.k2_for.sp12 = 0.95
species.bioen.forage.k2_for.sp13 = 0.95
species.bioen.forage.k2_for.sp14 = 0.95
species.bioen.forage.k2_for.sp15 = 1.1

############### alpha values (no unity)
species.beta.sp0 = 0.75
species.beta.sp1 = 0.75
species.beta.sp2 = 0.75
species.beta.sp3 = 0.75
species.beta.sp4 = 0.75
species.beta.sp5 = 0.75
species.beta.sp6 = 0.75
species.beta.sp7 = 0.75
species.beta.sp8 = 0.75
species.beta.sp9 = 0.75
species.beta.sp10 = 0.75
species.beta.sp11 = 0.75
species.beta.sp12 = 0.75
species.beta.sp13 = 0.75
species.beta.sp14 = 0.75
species.beta.sp15 = 0.75

############## r for kappa calculation (no unity)
# barrier.n: values set but fully checked.
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.r.sp0 into species.maturity.r.sp0
species.maturity.r.sp0 = 0.773584915
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.r.sp1 into species.maturity.r.sp1
species.maturity.r.sp1 = 0.747082654
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.r.sp2 into species.maturity.r.sp2
species.maturity.r.sp2 = 1.198753295
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.r.sp3 into species.maturity.r.sp3
species.maturity.r.sp3 = 0.558968521
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.r.sp4 into species.maturity.r.sp4
species.maturity.r.sp4 = 0.502980957
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.r.sp5 into species.maturity.r.sp5
species.maturity.r.sp5 = 0.519826119
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.r.sp6 into species.maturity.r.sp6
species.maturity.r.sp6 = 0.604330004
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.r.sp7 into species.maturity.r.sp7
species.maturity.r.sp7 = 0.208204975
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.r.sp8 into species.maturity.r.sp8
species.maturity.r.sp8 = 0.540768596
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.r.sp9 into species.maturity.r.sp9
species.maturity.r.sp9 = 0.585299763
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.r.sp10 into species.maturity.r.sp10
species.maturity.r.sp10 = 0.55339905
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.r.sp11 into species.maturity.r.sp11
species.maturity.r.sp11 = 0.728699757
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.r.sp12 into species.maturity.r.sp12
species.maturity.r.sp12 = 0.504756271
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.r.sp13 into species.maturity.r.sp13
species.maturity.r.sp13 = 0.866257872
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.r.sp14 into species.maturity.r.sp14
species.maturity.r.sp14 = 0.28643527
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.r.sp15 into species.maturity.r.sp15
species.maturity.r.sp15 = 1.65

species.maturity.eta.sp0 = 1
species.maturity.eta.sp1 = 1
species.maturity.eta.sp2 = 1
species.maturity.eta.sp3 = 1
species.maturity.eta.sp4 = 1
species.maturity.eta.sp5 = 1
species.maturity.eta.sp6 = 1
species.maturity.eta.sp7 = 1
species.maturity.eta.sp8 = 1
species.maturity.eta.sp9 = 1
species.maturity.eta.sp10 = 1
species.maturity.eta.sp11 = 1
species.maturity.eta.sp12 = 1
species.maturity.eta.sp13 = 1
species.maturity.eta.sp14 = 1
species.maturity.eta.sp15 = 1

### Maturité à Lm sans variabilité de Lm
# Slope (cm.L^-1)
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m0.sp0 into species.maturity.m0.sp0
species.maturity.m0.sp0 = 25.35853356
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m0.sp1 into species.maturity.m0.sp1
species.maturity.m0.sp1 = 25.76385255
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m0.sp2 into species.maturity.m0.sp2
species.maturity.m0.sp2 = 15
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m0.sp3 into species.maturity.m0.sp3
species.maturity.m0.sp3 = 12.95650186
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m0.sp4 into species.maturity.m0.sp4
species.maturity.m0.sp4 = 19.21510567
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m0.sp5 into species.maturity.m0.sp5
species.maturity.m0.sp5 = 43.513276
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m0.sp6 into species.maturity.m0.sp6
species.maturity.m0.sp6 = 19.170139
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m0.sp7 into species.maturity.m0.sp7
species.maturity.m0.sp7 = 74.09623922
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m0.sp8 into species.maturity.m0.sp8
species.maturity.m0.sp8 = 75.2574941
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m0.sp9 into species.maturity.m0.sp9
species.maturity.m0.sp9 = 32.85116924
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m0.sp10 into species.maturity.m0.sp10
species.maturity.m0.sp10 = 18.5
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m0.sp11 into species.maturity.m0.sp11
species.maturity.m0.sp11 = 21.72357681
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m0.sp12 into species.maturity.m0.sp12
species.maturity.m0.sp12 = 19.02561593
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m0.sp13 into species.maturity.m0.sp13
species.maturity.m0.sp13 = 28.89412443
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m0.sp14 into species.maturity.m0.sp14
species.maturity.m0.sp14 = 34.14893257
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m0.sp15 into species.maturity.m0.sp15
species.maturity.m0.sp15 = 5.37

# Intercepte (cm)
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m1.sp0 into species.maturity.m1.sp0
species.maturity.m1.sp0 = -0.89652575
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m1.sp1 into species.maturity.m1.sp1
species.maturity.m1.sp1 = 0.782315901
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m1.sp2 into species.maturity.m1.sp2
species.maturity.m1.sp2 = 0
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m1.sp3 into species.maturity.m1.sp3
species.maturity.m1.sp3 = 0
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m1.sp4 into species.maturity.m1.sp4
species.maturity.m1.sp4 = -3.544247862
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m1.sp5 into species.maturity.m1.sp5
species.maturity.m1.sp5 = -4.988393115
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m1.sp6 into species.maturity.m1.sp6
species.maturity.m1.sp6 = 0.693156414
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m1.sp7 into species.maturity.m1.sp7
species.maturity.m1.sp7 = -1.75140079
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m1.sp8 into species.maturity.m1.sp8
species.maturity.m1.sp8 = -2.744194529
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m1.sp9 into species.maturity.m1.sp9
species.maturity.m1.sp9 = -1.622630619
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m1.sp10 into species.maturity.m1.sp10
species.maturity.m1.sp10 = 0
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m1.sp11 into species.maturity.m1.sp11
species.maturity.m1.sp11 = -0.332275169
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m1.sp12 into species.maturity.m1.sp12
species.maturity.m1.sp12 = -2.474921686
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m1.sp13 into species.maturity.m1.sp13
species.maturity.m1.sp13 = -6.813415518
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m1.sp14 into species.maturity.m1.sp14
species.maturity.m1.sp14 = 1.472863302
# Osmose 4.4.1 - Renamed parameter species.bioen.maturity.m1.sp15 into species.maturity.m1.sp15
species.maturity.m1.sp15 = 0

predation.c.bioen.sp0 = 3.298118085
predation.c.bioen.sp1 = 3.710023779
predation.c.bioen.sp2 = 2.86
predation.c.bioen.sp3 = 1.509891027
predation.c.bioen.sp4 = 1.96
predation.c.bioen.sp5 = 2.690287905
predation.c.bioen.sp6 = 3.223440977
predation.c.bioen.sp7 = 3.408182295
predation.c.bioen.sp8 = 7.398372882
predation.c.bioen.sp9 = 4.642441914
predation.c.bioen.sp10 = 2.898103795
predation.c.bioen.sp11 = 4.732915648
predation.c.bioen.sp12 = 2.714847803
predation.c.bioen.sp13 = 3.639344385
predation.c.bioen.sp14 = 4.645947628
predation.c.bioen.sp15 = 2.28

species.larvae.growth.threshold.age.sp0 = 1
species.larvae.growth.threshold.age.sp1 = 1
species.larvae.growth.threshold.age.sp2 = 1
species.larvae.growth.threshold.age.sp3 = 1
species.larvae.growth.threshold.age.sp4 = 1
species.larvae.growth.threshold.age.sp5 = 1
species.larvae.growth.threshold.age.sp6 = 1
species.larvae.growth.threshold.age.sp7 = 1
species.larvae.growth.threshold.age.sp8 = 1
species.larvae.growth.threshold.age.sp9 = 1
species.larvae.growth.threshold.age.sp10 = 1
species.larvae.growth.threshold.age.sp11 = 1
species.larvae.growth.threshold.age.sp12 = 1
species.larvae.growth.threshold.age.sp13 = 1
species.larvae.growth.threshold.age.sp14 = 1
species.larvae.growth.threshold.age.sp15 = 1

predation.larval.ingestion.rate.increase.ratio.sp0 = 3.360509
predation.larval.ingestion.rate.increase.ratio.sp1 = 3.928367
predation.larval.ingestion.rate.increase.ratio.sp2 = 3.2
predation.larval.ingestion.rate.increase.ratio.sp3 = 5.524804
predation.larval.ingestion.rate.increase.ratio.sp4 = 5.664828
predation.larval.ingestion.rate.increase.ratio.sp5 = 4.463671
predation.larval.ingestion.rate.increase.ratio.sp6 = 2.852919
predation.larval.ingestion.rate.increase.ratio.sp7 = 4.470257
predation.larval.ingestion.rate.increase.ratio.sp8 = 2.102077
predation.larval.ingestion.rate.increase.ratio.sp9 = 2.914057
predation.larval.ingestion.rate.increase.ratio.sp10 = 3.360509
predation.larval.ingestion.rate.increase.ratio.sp11 = 2.480999
predation.larval.ingestion.rate.increase.ratio.sp12 = 3.201828
predation.larval.ingestion.rate.increase.ratio.sp13 = 2.789831
predation.larval.ingestion.rate.increase.ratio.sp14 = 4.114943
predation.larval.ingestion.rate.increase.ratio.sp15 = 1.282779

species.bioen.assimilation.sp0 = 0.8
species.bioen.assimilation.sp1 = 0.8
species.bioen.assimilation.sp2 = 0.8
species.bioen.assimilation.sp3 = 0.8
species.bioen.assimilation.sp4 = 0.8
species.bioen.assimilation.sp5 = 0.8
species.bioen.assimilation.sp6 = 0.8
species.bioen.assimilation.sp7 = 0.8
species.bioen.assimilation.sp8 = 0.8
species.bioen.assimilation.sp9 = 0.8
species.bioen.assimilation.sp10 = 0.8
species.bioen.assimilation.sp11 = 0.8
species.bioen.assimilation.sp12 = 0.8
species.bioen.assimilation.sp13 = 0.8
species.bioen.assimilation.sp14 = 0.8
species.bioen.assimilation.sp15 = 0.8

######################### physical dataset parameters
temperature.varname = T
temperature.filename = physical_ns.nc
temperature.factor = 1
temperature.offset = 0
temperature.nsteps.year = 24

oxygen.varname = O2
oxygen.filename = physical_ns.nc
# conversion of O2 from mmol/m3 into mmol/L,
# conversion of mmol/m3 into mg/L is hard coded in Osmose,
# conversion mmol/m3 -> mmol/L: *= 1e-3,
# O2 en %saturation --> pas besoin de factor ,
oxygen.factor = 1
oxygen.offset = 0
oxygen.nsteps.year = 24

# Trait parameters -------------------------------------------------------------

# Osmose 4.4.0 - Renamed parameter simulation.genetic.enabled into module.genetics.enabled
module.genetics.enabled = TRUE
output.evolvingtraits.enabled = TRUE
population.genotype.transmission.year.start = 50

output.allele.frequency.enabled =TRUE
output.expected.heterozygosity.enabled = TRUE
output.observed.heterozygosity.enabled = TRUE
output.observed.heterozygosity.byage.enabled = TRUE
output.allele.frequency.byage.enabled = TRUE
output.expected.heterozygosity.byage.enabled = TRUE

neutral.trait.mean.sp0=0
neutral.trait.mean.sp1=0
neutral.trait.mean.sp2=0
neutral.trait.mean.sp3=0
neutral.trait.mean.sp4=0
neutral.trait.mean.sp5=0
neutral.trait.mean.sp6=0
neutral.trait.mean.sp7=0
neutral.trait.mean.sp8=0
neutral.trait.mean.sp9=0
neutral.trait.mean.sp10=0
neutral.trait.mean.sp11=0
neutral.trait.mean.sp12=0
neutral.trait.mean.sp13=0
neutral.trait.mean.sp14=0
neutral.trait.mean.sp15=0

neutral.trait.nlocus.sp0=100
neutral.trait.nlocus.sp1=100
neutral.trait.nlocus.sp2=100
neutral.trait.nlocus.sp3=100
neutral.trait.nlocus.sp4=100
neutral.trait.nlocus.sp5=100
neutral.trait.nlocus.sp6=100
neutral.trait.nlocus.sp7=100
neutral.trait.nlocus.sp8=100
neutral.trait.nlocus.sp9=100
neutral.trait.nlocus.sp10=100
neutral.trait.nlocus.sp11=100
neutral.trait.nlocus.sp12=100
neutral.trait.nlocus.sp13=100
neutral.trait.nlocus.sp14=100
neutral.trait.nlocus.sp15=100

neutral.trait.nval.sp0=2
neutral.trait.nval.sp1=2
neutral.trait.nval.sp2=2
neutral.trait.nval.sp3=2
neutral.trait.nval.sp4=2
neutral.trait.nval.sp5=2
neutral.trait.nval.sp6=2
neutral.trait.nval.sp7=2
neutral.trait.nval.sp8=2
neutral.trait.nval.sp9=2
neutral.trait.nval.sp10=2
neutral.trait.nval.sp11=2
neutral.trait.nval.sp12=2
neutral.trait.nval.sp13=2
neutral.trait.nval.sp14=2
neutral.trait.nval.sp15=2

neutral.trait.var.sp0=0.2
neutral.trait.var.sp1=0.2
neutral.trait.var.sp2=0.2
neutral.trait.var.sp3=0.2
neutral.trait.var.sp4=0.2
neutral.trait.var.sp5=0.2
neutral.trait.var.sp6=0.2
neutral.trait.var.sp7=0.2
neutral.trait.var.sp8=0.2
neutral.trait.var.sp9=0.2
neutral.trait.var.sp10=0.2
neutral.trait.var.sp11=0.2
neutral.trait.var.sp12=0.2
neutral.trait.var.sp13=0.2
neutral.trait.var.sp14=0.2
neutral.trait.var.sp15=0.2

neutral.trait.envvar.sp0=0
neutral.trait.envvar.sp1=0
neutral.trait.envvar.sp2=0
neutral.trait.envvar.sp3=0
neutral.trait.envvar.sp4=0
neutral.trait.envvar.sp5=0
neutral.trait.envvar.sp6=0
neutral.trait.envvar.sp7=0
neutral.trait.envvar.sp8=0
neutral.trait.envvar.sp9=0
neutral.trait.envvar.sp10=0
neutral.trait.envvar.sp11=0
neutral.trait.envvar.sp12=0
neutral.trait.envvar.sp13=0
neutral.trait.envvar.sp14=0
neutral.trait.envvar.sp15=0
# 
# imax.trait.mean.sp0 = 14.0161771995127
# imax.trait.mean.sp1 = 16.4430568215567
# imax.trait.mean.sp2 = 9.48185907748328
# imax.trait.mean.sp3 = 12.119810661574
# imax.trait.mean.sp4 = 9.91981236341634
# imax.trait.mean.sp5 = 10.9236872467061
# imax.trait.mean.sp6 = 9.94816846789243
# imax.trait.mean.sp7 = 14.4507137847945
# imax.trait.mean.sp8 = 21.5633912243199
# imax.trait.mean.sp9 = 15.7906712579798
# imax.trait.mean.sp10 = 13.3885972672411
# imax.trait.mean.sp11 = 17.8426721069184
# imax.trait.mean.sp12 = 8.78878491915743
# imax.trait.mean.sp13 = 13.5696659330434
# imax.trait.mean.sp14 = 17.2822647633446
# imax.trait.mean.sp15 = 8.91297095163207
# 
# #var.genet = (0.01*Imax)²
# imax.trait.var.sp0 = 0
# imax.trait.var.sp1 = 0
# imax.trait.var.sp2 = 0
# imax.trait.var.sp3 = 0
# imax.trait.var.sp4 = 0
# imax.trait.var.sp5 = 0
# imax.trait.var.sp6 = 0
# imax.trait.var.sp7 = 0
# imax.trait.var.sp8 = 0
# imax.trait.var.sp9 = 0
# imax.trait.var.sp10 = 0
# imax.trait.var.sp11 = 0
# imax.trait.var.sp12 = 0
# imax.trait.var.sp13 = 0
# imax.trait.var.sp14 = 0
# imax.trait.var.sp15 = 0
# 
# imax.trait.nlocus.sp0 = 10
# imax.trait.nlocus.sp1 = 10
# imax.trait.nlocus.sp2 = 10
# imax.trait.nlocus.sp3 = 10
# imax.trait.nlocus.sp4 = 10
# imax.trait.nlocus.sp5 = 10
# imax.trait.nlocus.sp6 = 10
# imax.trait.nlocus.sp7 = 10
# imax.trait.nlocus.sp8 = 10
# imax.trait.nlocus.sp9 = 10
# imax.trait.nlocus.sp10 = 10
# imax.trait.nlocus.sp11 = 10
# imax.trait.nlocus.sp12 = 10
# imax.trait.nlocus.sp13 = 10
# imax.trait.nlocus.sp14 = 10
# imax.trait.nlocus.sp15 = 10
# 
# imax.trait.nval.sp0 = 7
# imax.trait.nval.sp1 = 7
# imax.trait.nval.sp2 = 7
# imax.trait.nval.sp3 = 7
# imax.trait.nval.sp4 = 7
# imax.trait.nval.sp5 = 7
# imax.trait.nval.sp6 = 7
# imax.trait.nval.sp7 = 7
# imax.trait.nval.sp8 = 7
# imax.trait.nval.sp9 = 7
# imax.trait.nval.sp10 = 7
# imax.trait.nval.sp11 = 7
# imax.trait.nval.sp12 = 7
# imax.trait.nval.sp13 = 7
# imax.trait.nval.sp14 = 7
# imax.trait.nval.sp15 = 7
# 
# imax.trait.envvar.sp0 = 0
# imax.trait.envvar.sp1 = 0
# imax.trait.envvar.sp2 = 0
# imax.trait.envvar.sp3 = 0
# imax.trait.envvar.sp4 = 0
# imax.trait.envvar.sp5 = 0
# imax.trait.envvar.sp6 = 0
# imax.trait.envvar.sp7 = 0
# imax.trait.envvar.sp8 = 0
# imax.trait.envvar.sp9 = 0
# imax.trait.envvar.sp10 = 0
# imax.trait.envvar.sp11 = 0
# imax.trait.envvar.sp12 = 0
# imax.trait.envvar.sp13 = 0
# imax.trait.envvar.sp14 = 0
# imax.trait.envvar.sp15 = 0

# Fishery parameters -----------------------------------------------------------

# Osmose 4.4.0 -   Did not rename fisheries.enabled into process.multispecies.fisheries.enabled as it is already defined in fileR:\sync\Github\_temp\osmose-ns\osmose-ns-v4.4.1_clim\New_param_Eco_Bioen.R
# fisheries.enabled = FALSE
# Osmose 4.4.0 - Renamed parameter process.multispecies.fisheries.enabled into module.multispecies.fisheries.enabled
module.multispecies.fisheries.enabled = FALSE
#fisheries.check.enabled = TRUE
#simulation.fishing.mortality.enabled = TRUE
#simulation.nfisheries = 16
#fisheries.movement.netcdf.enabled = FALSE

#fisheries.catchability.file = Fishery_input/ns_catchability_sp.csv
#fisheries.discards.file = Fishery_input/ns_discards_sp.csv

# Fisheries = 1 per species
#fisheries.name.fsh0 = FAtlanticHerring
#fisheries.name.fsh1 = FAtlanticMackerel
#fisheries.name.fsh2 = FSandeel
#fisheries.name.fsh3 = FEuropeanSprat
#fisheries.name.fsh4 = FNorwayPout
#fisheries.name.fsh5 = FEuropeanPlaice
#fisheries.name.fsh6 = FCommonSole
#fisheries.name.fsh7 = FSaithe
#fisheries.name.fsh8 = FAtlanticCod
#fisheries.name.fsh9 = FHaddock
#fisheries.name.fsh10 = FHorseMackerel
#fisheries.name.fsh11 = FWhiting
#fisheries.name.fsh12 = FCommonDab
#fisheries.name.fsh13 = FGreyGurnard
#fisheries.name.fsh14 = FHake
#fisheries.name.fsh15 = FShrimp

# Economy parameters -----------------------------------------------------------

# Osmose 4.4.0 - Renamed parameter economy.enabled into module.bioeconomics.enabled
module.bioeconomics.enabled = FALSE

economic.output.stage.structure = weight

economic.output.stage.threshold.sp0 = 31, 57, 85, 125, 250
economic.output.stage.threshold.sp1 = 80, 100, 200, 500
economic.output.stage.threshold.sp2 = 1
economic.output.stage.threshold.sp3 = 4
economic.output.stage.threshold.sp4 = 5
economic.output.stage.threshold.sp5 = 150, 300, 400, 600
economic.output.stage.threshold.sp6 = 120, 170, 250, 330, 500
economic.output.stage.threshold.sp7 = 300, 1500, 3000, 5000
economic.output.stage.threshold.sp8 = 300, 1000, 2000, 4000, 7000
economic.output.stage.threshold.sp9 = 170, 370, 570, 1000
economic.output.stage.threshold.sp10 = 20, 80, 200, 400, 600
economic.output.stage.threshold.sp11 = 110, 250, 350, 500
economic.output.stage.threshold.sp12 = 130, 250
economic.output.stage.threshold.sp13 = 60, 200
economic.output.stage.threshold.sp14 = 150, 280, 600, 1200, 2500
economic.output.stage.threshold.sp15 = 1

#elasticity of substitution between species
substitution.elasticity = 4

### stock elasticity <1 means stock hyperstable -> fix to 0.8 for now
species.stock.elasticity.sp0 = 0.0
species.stock.elasticity.sp1 = 2.460245e-01
species.stock.elasticity.sp2 = 0.0
species.stock.elasticity.sp3 = 2.460665e-01
species.stock.elasticity.sp4 = 0.010796995
species.stock.elasticity.sp5 = 0.0
species.stock.elasticity.sp6 = 6.285602e-02
species.stock.elasticity.sp7 = 0.01119766
species.stock.elasticity.sp8 = 1.582588e-01
species.stock.elasticity.sp9 = 0.0
species.stock.elasticity.sp10 = 0.0
species.stock.elasticity.sp11 = 2.046791e-01
species.stock.elasticity.sp12 = 0.0
species.stock.elasticity.sp13 = 0.04229606
species.stock.elasticity.sp14 = 0.0
species.stock.elasticity.sp15 = 7.325882e-02

# Baseline costs at t0
baseline.costs.t0.sp0 = 326.4156126
baseline.costs.t0.sp1 = 1.728810e+04
baseline.costs.t0.sp2 = 2.972616e+02
baseline.costs.t0.sp3 = 8.149485e+03
baseline.costs.t0.sp4 = 140.777584402
baseline.costs.t0.sp5 = 1022.30977539
baseline.costs.t0.sp6 = 1.459262e+04
baseline.costs.t0.sp7 = 1100.90821928
baseline.costs.t0.sp8 = 1.460089e+04
baseline.costs.t0.sp9 = 838.62931302
baseline.costs.t0.sp10 = 319.02630279
baseline.costs.t0.sp11 = 9.229539e+03
baseline.costs.t0.sp12 = 586.273208434
baseline.costs.t0.sp13 = 969.58568049
baseline.costs.t0.sp14 = 1.223151e+03
baseline.costs.t0.sp15 = 5.541426e+03

size.consumption.elasticity = 4

species.sizeconsumption.elasticity.sp0 = 8
species.sizeconsumption.elasticity.sp1 = 8
species.sizeconsumption.elasticity.sp2 = 3
species.sizeconsumption.elasticity.sp3 = 3
species.sizeconsumption.elasticity.sp4 = 3
species.sizeconsumption.elasticity.sp5 = 8
species.sizeconsumption.elasticity.sp6 = 8
species.sizeconsumption.elasticity.sp7 = 8
species.sizeconsumption.elasticity.sp8 = 8
species.sizeconsumption.elasticity.sp9 = 8
species.sizeconsumption.elasticity.sp10 = 8
species.sizeconsumption.elasticity.sp11 = 8
species.sizeconsumption.elasticity.sp12 = 8
species.sizeconsumption.elasticity.sp13 = 8
species.sizeconsumption.elasticity.sp14 = 8
species.sizeconsumption.elasticity.sp15 = 3

species.consumption.elasticity.sp0 = 0.045215600
species.consumption.elasticity.sp1 = 0.052296754
species.consumption.elasticity.sp2 = 0.004765379
species.consumption.elasticity.sp3 = 0.004367975
species.consumption.elasticity.sp4 = 0.001166938
species.consumption.elasticity.sp5 = 0.051370748
species.consumption.elasticity.sp6 = 0.429237064
species.consumption.elasticity.sp7 = 0.042518198
species.consumption.elasticity.sp8 = 0.163673297
species.consumption.elasticity.sp9 = 0.035624833
species.consumption.elasticity.sp10 = 0.019957743
species.consumption.elasticity.sp11 = 0.043043137
species.consumption.elasticity.sp12 = 0.003428198
species.consumption.elasticity.sp13 = 0.002608459
species.consumption.elasticity.sp14 = 0.086545029
species.consumption.elasticity.sp15 = 0.014180650

price.time.trend.sp0 = 0.0147805
price.time.trend.sp1 = 3.576891e-02
price.time.trend.sp2 = 1.896321e-04
price.time.trend.sp3 = -8.880207e-03
price.time.trend.sp4 = -0.009320397
price.time.trend.sp5 = -0.02262025
price.time.trend.sp6 = -1.723837e-02
price.time.trend.sp7 = -0.01045096
price.time.trend.sp8 = -2.864259e-03
price.time.trend.sp9 = -0.00458941
price.time.trend.sp10 = 0.01828865
price.time.trend.sp11 = 6.768317e-04
price.time.trend.sp12 = -0.008905486
price.time.trend.sp13 = -0.03155171
price.time.trend.sp14 = 5.755827e-03
price.time.trend.sp15 = 1.554470e-02

species.size.preference.sp0 = 0.1436690, 0.2432393, 0.2429449, 0.1940378, 0.1761089
species.size.preference.sp1 = 0.23715467, 0.26470177, 0.41424237, 0.08390119
species.size.preference.sp2 = 1
species.size.preference.sp3 = 1
species.size.preference.sp4 = 1
species.size.preference.sp5 = 0.2089878, 0.2698584, 0.2653381, 0.2558157
species.size.preference.sp6 = 0.1929576, 0.1968394, 0.2124119, 0.2269605, 0.1708306
species.size.preference.sp7 = 0.2692069, 0.4055411, 0.2094138, 0.1158381
species.size.preference.sp8 = 0.1192954, 0.1854998, 0.2561656, 0.2318471, 0.2071920
species.size.preference.sp9 = 0.2647919, 0.2856680, 0.3027170, 0.1468232
species.size.preference.sp10 = 0.19446340, 0.26600924, 0.27504337, 0.16471691, 0.09976708
species.size.preference.sp11 = 0.2081937, 0.2263748, 0.2969661, 0.2684654
species.size.preference.sp12 = 0.4146418, 0.5853582
species.size.preference.sp13 = 0.7878125, 0.2121875
species.size.preference.sp14 = 0.1441224, 0.1616195, 0.2417079, 0.2461011, 0.2064491
species.size.preference.sp15 = 1

elasticity.demand.fish = 4.480645

weight.fish.consumption = 552723.2
