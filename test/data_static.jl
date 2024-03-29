test_LP() = CoreModel(
    zeros(4, 3),
    zeros(4),
    ones(3),
    ones(3),
    ones(3),
    ["r$x" for x = 1:3],
    ["m$x" for x = 1:4],
)

test_simpleLP() = CoreModel(
    [
        1.0 1.0
        -1.0 1.0
    ],
    [3.0, 1.0],
    [-0.25, 1.0],
    -ones(2),
    2.0 * ones(2),
    ["r$x" for x = 1:2],
    ["m$x" for x = 1:2],
)

test_simpleLP2() = CoreModel(
    zeros(2, 2),
    [0.0, 0.0],
    [-0.25, 1.0],
    -ones(2),
    2.0 * ones(2),
    ["r$x" for x = 1:2],
    ["m$x" for x = 1:2],
)

test_sparseLP() = CoreModel(
    sprand(4000, 3000, 0.5),
    sprand(4000, 0.5),
    sprand(3000, 0.5),
    sprand(3000, 0.5),
    sprand(3000, 0.5),
    ["r$x" for x = 1:3000],
    ["m$x" for x = 1:4000],
)

test_coupledLP() = CoreModelCoupled(
    CoreModel(
        sprand(4000, 3000, 0.5),
        sprand(4000, 0.5),
        sprand(3000, 0.5),
        sprand(3000, 0.5),
        sprand(3000, 0.5),
        ["r$x" for x = 1:3000],
        ["m$x" for x = 1:4000],
    ),
    sprand(2000, 3000, 0.5),
    sprand(2000, 0.5),
    sprand(2000, 0.5),
)

test_toyModel() = CoreModel(
    [
        -1.0 1.0 0.0 0.0 0.0 0.0 0.0
        -2.0 0.0 1.0 0.0 0.0 0.0 0.0
        1.0 0.0 0.0 0.0 0.0 0.0 -1.0
        0.0 -1.0 0.0 -1.0 0.0 0.0 0.0
        0.0 0.0 -1.0 0.0 -1.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 -1.0 1.0
    ],
    zeros(6),
    [0, 0, 0, 0, 0, 0, 1.0],
    fill(-1000.0, 7),
    fill(1000.0, 7),
    ["r1", "m1t", "m3t", "EX_m1(e)", "EX_m3(e)", "EX_biomass(e)", "biomass1"],
    ["m1[c]", "m3[c]", "m2[c]", "m1[e]", "m3[e]", "biomass[c]"],
)

const reaction_standard_gibbs_free_energies = Dict{String,Float64}(
    #=
    ΔᵣG⁰ data from Equilibrator using the E. coli core model's reactions
    To generate this data manually, go to https://equilibrator.weizmann.ac.il/ and
    enter each reaction into the search bar, fill in the ΔG⁰ below from there. To generate
    automatically, use the eQuilibrator.jl package.
    =#
    "ACALD" => -21.268198981756314,
    "PTAr" => 8.651025243027988,
    "ALCD2x" => 17.47260052762408,
    "PDH" => -34.246780702043225,
    "PYK" => -24.48733600711958,
    "CO2t" => 0.0,
    "MALt2_2" => -6.839209921974724,
    "CS" => -39.330940148207475,
    "PGM" => -4.470553692565886,
    "TKT1" => -1.4976299544215408,
    "ACONTa" => 8.46962176350985,
    "GLNS" => -15.771242654033706,
    "ICL" => 9.53738025507404,
    "FBA" => 23.376920310319235,
    "SUCCt3" => -43.97082007726285,
    "FORt2" => -3.4211767267069124,
    "G6PDH2r" => -7.3971867270165035,
    "AKGDH" => -28.235442362320782,
    "TKT2" => -10.318436817276165,
    "FRD7" => 73.61589524884772,
    "SUCOAS" => -1.1586968559555615,
    "FBP" => -11.606887851480572,
    "ICDHyr" => 5.398885342794983,
    "AKGt2r" => 10.085299238039797,
    "GLUSy" => -47.21690395283849,
    "TPI" => 5.621932460512994,
    "FORt" => 13.509690780237293,
    "ACONTb" => -1.622946931741609,
    "GLNabc" => -30.194559792842753,
    "RPE" => -3.3881719029747615,
    "ACKr" => 14.027197131450492,
    "THD2" => -33.84686533309243,
    "PFL" => -19.814421615735,
    "RPI" => 4.477649590945703,
    "D_LACt2" => -3.4223903975852004,
    "TALA" => -0.949571985515206,
    "PPCK" => 10.659841402564751,
    "ACt2r" => -3.4196348363397995,
    "NH4t" => -13.606633927039097,
    "PGL" => -25.94931748161696,
    "NADTRHD" => -0.014869680795754903,
    "PGK" => 19.57192102020454,
    "LDH_D" => 20.04059765689044,
    "ME1" => 12.084968268076864,
    "PIt2r" => 10.415108493818785,
    "ATPS4r" => -37.570267233299816,
    "PYRt2" => -3.422891289768689,
    "GLCpts" => -45.42430981510088,
    "GLUDy" => 32.834943812395665,
    "CYTBD" => -59.700410815493775,
    "FUMt2_2" => -6.845105500839001,
    "FRUpts2" => -42.67529760694199,
    "GAPD" => 0.5307809794271634,
    "H2Ot" => -1.5987211554602254e-14,
    "PPC" => -40.81304419704113,
    "NADH16" => -80.37770501380615,
    "PFK" => -18.546314942995934,
    "MDH" => 25.912462872631522,
    "PGI" => 2.6307087407442395,
    "O2t" => 0.0,
    "ME2" => 12.099837948872533,
    "GND" => 10.312275879236381,
    "SUCCt2_2" => -6.82178244356977,
    "GLUN" => -14.381960140443113,
    "ETOHt2r" => -16.930867506944217,
    "ADK1" => 0.3893321583896068,
    "ACALDt" => -3.197442310920451e-14,
    "SUCDi" => -73.61589524884772,
    "ENO" => -3.8108376097261782,
    "MALS" => -39.22150045995042,
    "GLUt2r" => -3.499043558772904,
    "PPS" => -6.0551989457468665,
    "FUM" => -3.424133018702122,
)

const ecoli_core_gene_product_masses = Dict{String,Float64}(
    #=
    Data downloaded from Uniprot for E. coli K12,
    gene mass in kDa. To obtain these data yourself, go to
    Uniprot: https://www.uniprot.org/
    and search using these terms: <reviewed:yes AND organism:"Escherichia coli (strain K12) [83333]">
    =#
    "b4301" => 23.214,
    "b1602" => 48.723,
    "b4154" => 65.972,
    "b3236" => 32.337,
    "b1621" => 56.627,
    "b1779" => 35.532,
    "b3951" => 85.96,
    "b1676" => 50.729,
    "b3114" => 85.936,
    "b1241" => 96.127,
    "b2276" => 52.044,
    "b1761" => 48.581,
    "b3925" => 35.852,
    "b3493" => 53.389,
    "b3733" => 31.577,
    "b2926" => 41.118,
    "b0979" => 42.424,
    "b4015" => 47.522,
    "b2296" => 43.29,
    "b4232" => 36.834,
    "b3732" => 50.325,
    "b2282" => 36.219,
    "b2283" => 100.299,
    "b0451" => 44.515,
    "b2463" => 82.417,
    "b0734" => 42.453,
    "b3738" => 30.303,
    "b3386" => 24.554,
    "b3603" => 59.168,
    "b2416" => 63.562,
    "b0729" => 29.777,
    "b0767" => 36.308,
    "b3734" => 55.222,
    "b4122" => 60.105,
    "b2987" => 53.809,
    "b2579" => 14.284,
    "b0809" => 26.731,
    "b1524" => 33.516,
    "b3612" => 56.194,
    "b3735" => 19.332,
    "b3731" => 15.068,
    "b1817" => 35.048,
    "b1603" => 54.623,
    "b1773" => 30.81,
    "b4090" => 16.073,
    "b0114" => 99.668,
    "b3962" => 51.56,
    "b2464" => 35.659,
    "b2976" => 80.489,
    "b1818" => 27.636,
    "b2285" => 18.59,
    "b1702" => 87.435,
    "b1849" => 42.434,
    "b1812" => 50.97,
    "b0902" => 28.204,
    "b3403" => 59.643,
    "b1612" => 60.299,
    "b1854" => 51.357,
    "b0811" => 27.19,
    "b0721" => 14.299,
    "b2914" => 22.86,
    "b1297" => 53.177,
    "b0723" => 64.422,
    "b3919" => 26.972,
    "b3115" => 43.384,
    "b4077" => 47.159,
    "b3528" => 45.436,
    "b0351" => 33.442,
    "b2029" => 51.481,
    "b1819" => 30.955,
    "b0728" => 41.393,
    "b2935" => 72.212,
    "b2415" => 9.119,
    "b0727" => 44.011,
    "b0116" => 50.688,
    "b0485" => 32.903,
    "b3736" => 17.264,
    "b0008" => 35.219,
    "b3212" => 163.297,
    "b3870" => 51.904,
    "b4014" => 60.274,
    "b2280" => 19.875,
    "b2133" => 64.612,
    "b2278" => 66.438,
    "b0118" => 93.498,
    "b2288" => 16.457,
    "b3739" => 13.632,
    "b3916" => 34.842,
    "b3952" => 32.43,
    "b2925" => 39.147,
    "b2465" => 73.043,
    "b2297" => 77.172,
    "b2417" => 18.251,
    "b4395" => 24.065,
    "b3956" => 99.063,
    "b0722" => 12.868,
    "b2779" => 45.655,
    "b0115" => 66.096,
    "b0733" => 58.205,
    "b1478" => 35.38,
    "b2492" => 30.565,
    "b0724" => 26.77,
    "b0755" => 28.556,
    "b1136" => 45.757,
    "b2286" => 68.236,
    "b0978" => 57.92,
    "b1852" => 55.704,
    "b2281" => 20.538,
    "b2587" => 47.052,
    "b2458" => 36.067,
    "b0904" => 30.991,
    "b1101" => 50.677,
    "b0875" => 23.703,
    "b3213" => 52.015,
    "b2975" => 58.92,
    "b0720" => 48.015,
    "b0903" => 85.357,
    "b1723" => 32.456,
    "b2097" => 38.109,
    "b3737" => 8.256,
    "b0810" => 24.364,
    "b4025" => 61.53,
    "b1380" => 36.535,
    "b0356" => 39.359,
    "b2277" => 56.525,
    "b1276" => 97.677,
    "b4152" => 15.015,
    "b1479" => 63.197,
    "b4153" => 27.123,
    "b4151" => 13.107,
    "b2287" => 25.056,
    "b0474" => 23.586,
    "b2284" => 49.292,
    "b1611" => 50.489,
    "b0726" => 105.062,
    "b2279" => 10.845,
)

const ecoli_core_reaction_kcats = Dict{String,Vector{Tuple{Float64,Float64}}}(
    #=
    Data taken from Heckmann, David, et al. "Machine learning applied to enzyme
    turnover numbers reveals protein structural correlates and improves metabolic
    models." Nature communications 9.1 (2018): 1-10.
    =#
    "ACALD" =>
        [(568.1130792316333, 568.1130792316333), (568.856126503717, 568.856126503717)],
    "PTAr" => [
        (1171.9703624351055, 1171.9703624351055),
        (1173.7231032615289, 1173.7231032615289),
    ],
    "ALCD2x" => [
        (75.9547881894345, 75.9547881894345),
        (75.96334310351442, 75.96334310351442),
        (76.1472359297987, 76.1472359297987),
    ],
    "PDH" => [(529.7610874857239, 529.7610874857239)],
    "PYK" => [
        (422.0226052080562, 422.0226052080562),
        (422.1332899347833, 422.1332899347833),
    ],
    "MALt2_2" => [(234.03664660088714, 234.03664660088714)],
    "CS" => [(113.29607453875758, 113.29607453875758)],
    "PGM" => [
        (681.4234715886669, 681.4234715886669),
        (681.6540601244343, 681.6540601244343),
        (680.5234799168278, 680.5234799168278),
    ],
    "TKT1" => [
        (311.16139580671637, 311.16139580671637),
        (311.20967965149947, 311.20967965149947),
    ],
    "ACONTa" => [
        (191.02308213992006, 191.02308213992006),
        (191.03458045697235, 191.03458045697235),
    ],
    "GLNS" => [
        (89.83860937287024, 89.83860937287024),
        (89.82177852142014, 89.82177852142014),
    ],
    "ICL" => [(17.45922330097792, 17.45922330097792)],
    "FBA" => [
        (373.425646787578, 373.425646787578),
        (372.74936053215833, 372.74936053215833),
        (372.88627228768166, 372.88627228768166),
    ],
    "FORt2" => [
        (233.93045260179326, 233.93045260179326),
        (233.84804009142908, 233.84804009142908),
    ],
    "G6PDH2r" => [(589.3761070080022, 589.3761070080022)],
    "AKGDH" => [(264.48071159327156, 264.48071159327156)],
    "TKT2" => [
        (467.4226876901618, 467.4226876901618),
        (468.1440593542596, 468.1440593542596),
    ],
    "FRD7" => [(90.20637824912605, 90.20637824912605)],
    "SUCOAS" => [(18.494387648707622, 18.494387648707622)],
    "FBP" => [
        (568.5346256470805, 568.5346256470805),
        (567.6367759041788, 567.6367759041788),
    ],
    "ICDHyr" => [(39.62446791678959, 39.62446791678959)],
    "AKGt2r" => [(234.99097804446805, 234.99097804446805)],
    "GLUSy" => [(33.262997317319055, 33.262997317319055)],
    "TPI" => [(698.301904211076, 698.301904211076)],
    "FORt" => [
        (234.38391855848187, 234.38391855848187),
        (234.34725576182922, 234.34725576182922),
    ],
    "ACONTb" => [
        (159.74612206327865, 159.74612206327865),
        (159.81975755249232, 159.81975755249232),
    ],
    "GLNabc" => [(233.80358131677775, 233.80358131677775)],
    "RPE" => [
        (1772.4850826683305, 1772.4850826683305),
        (1768.8536177485582, 1768.8536177485582),
    ],
    "ACKr" => [
        (554.611547307207, 554.611547307207),
        (555.112707891257, 555.112707891257),
        (555.2464368932744, 555.2464368932744),
    ],
    "THD2" => [(24.739139801185537, 24.739139801185537)],
    "PFL" => [
        (96.56316095411077, 96.56316095411077),
        (96.65024313036014, 96.65024313036014),
        (96.60761818004025, 96.60761818004025),
        (96.49541118899961, 96.49541118899961),
    ],
    "RPI" => [
        (51.771578021074234, 51.771578021074234),
        (51.81603467243345, 51.81603467243345),
    ],
    "D_LACt2" => [
        (233.51709131524734, 233.51709131524734),
        (233.83187606098016, 233.83187606098016),
    ],
    "TALA" => [
        (109.05210545422884, 109.05210545422884),
        (109.04246437049026, 109.04246437049026),
    ],
    "PPCK" => [(218.4287805666016, 218.4287805666016)],
    "PGL" => [(2120.4297518987964, 2120.4297518987964)],
    "NADTRHD" => [
        (186.99387360624777, 186.99387360624777),
        (187.16629305266423, 187.16629305266423),
    ],
    "PGK" => [(57.641966636896335, 57.641966636896335)],
    "LDH_D" => [
        (31.11118891764946, 31.11118891764946),
        (31.12493425054357, 31.12493425054357),
    ],
    "ME1" => [(487.0161203971232, 487.0161203971232)],
    "PIt2r" => [
        (233.8651331835765, 233.8651331835765),
        (234.27374798581067, 234.27374798581067),
    ],
    "ATPS4r" => [
        (7120.878030435999, 7120.878030435999),
        (7116.751386037507, 7116.751386037507),
    ],
    "GLCpts" => [
        (233.9009878400008, 233.9009878400008),
        (233.66656882114864, 233.66656882114864),
        (233.66893882934883, 233.66893882934883),
    ],
    "GLUDy" => [(105.32811069172409, 105.32811069172409)],
    "CYTBD" => [
        (153.18512795009505, 153.18512795009505),
        (153.2429537682265, 153.2429537682265),
    ],
    "FUMt2_2" => [(234.37495609395967, 234.37495609395967)],
    "FRUpts2" => [(234.1933863380989, 234.1933863380989)],
    "GAPD" => [(128.76795529111456, 128.76795529111456)],
    "PPC" => [(165.52424516841342, 165.52424516841342)],
    "NADH16" => [(971.7487306963936, 971.7487306963936)],
    "PFK" => [
        (1000.4626204522712, 1000.4626204522712),
        (1000.5875517343595, 1000.5875517343595),
    ],
    "MDH" => [(25.931655783969283, 25.931655783969283)],
    "PGI" => [(468.11833198138834, 468.11833198138834)],
    "ME2" => [(443.0973626307168, 443.0973626307168)],
    "GND" => [(240.1252264230952, 240.1252264230952)],
    "SUCCt2_2" => [(234.18109388303225, 234.18109388303225)],
    "GLUN" => [
        (44.76358496525738, 44.76358496525738),
        (44.84850207360875, 44.84850207360875),
        (44.76185250415503, 44.76185250415503),
    ],
    "ADK1" => [(111.64869652600649, 111.64869652600649)],
    "SUCDi" => [(680.3193833053011, 680.3193833053011)],
    "ENO" => [(209.35855069219886, 209.35855069219886)],
    "MALS" => [
        (252.7540503869977, 252.7540503869977),
        (252.2359738678874, 252.2359738678874),
    ],
    "GLUt2r" => [(234.22890837451837, 234.22890837451837)],
    "PPS" => [(706.1455885214322, 706.1455885214322)],
    "FUM" => [
        (1576.8372583425075, 1576.8372583425075),
        (1576.233088455828, 1576.233088455828),
        (1575.9638204848736, 1575.9638204848736),
    ],
)
