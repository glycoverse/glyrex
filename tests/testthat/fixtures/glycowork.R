list(
  list(
    pattern = "Hex-HexNAc-([Hex|Fuc])?-HexNAc",
    glycan = "GalNAc(b1-4)GlcNAc(b1-2)Man(a1-3)[Neu5Gc(a2-6)GalNAc(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc",
    expected = list("Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc"),
    line = 3626L
  ),
  list(
    pattern = "Hex-HexNAc-([Hex|Fuc])*-HexNAc",
    glycan = "GalNAc(b1-4)GlcNAc(b1-2)Man(a1-3)[Neu5Gc(a2-6)GalNAc(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc",
    expected = list("Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc"),
    line = 3627L
  ),
  list(
    pattern = ".-.-([Hex|Fuc])+-.",
    glycan = "GalNAc(b1-4)GlcNAc(b1-2)Man(a1-3)[Neu5Gc(a2-6)GalNAc(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc",
    expected = list(
      "Neu5Gc(a2-6)GalNAc(b1-4)[Fuc(a1-3)]GlcNAc",
      "Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc",
      "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man",
      "GlcNAc(b1-2)Man(a1-6)[Man(a1-3)]Man"
    ),
    line = 3628L
  ),
  list(
    pattern = "Fuc-Galb3/4-([Hex|Fuc])*-HexNAc",
    glycan = "Fuc(a1-2)Gal(b1-?)[Fuc(a1-?)]GlcNAc(b1-6)[Gal(b1-3)]GalNAc",
    expected = list("Fuc(a1-2)Gal(b1-?)[Fuc(a1-?)]GlcNAc"),
    line = 3629L
  ),
  list(
    pattern = "Fuc-([^Gal])+-GlcNAc",
    glycan = "Fuc(a1-3)[GlcNAc(b1-4)]GlcNAc",
    expected = list(
      "Fuc(a1-3)[GlcNAc(b1-4)]GlcNAc"
    ),
    line = 3630L
  ),
  list(
    pattern = "Fuc-([^Gal])+-GlcNAc",
    glycan = "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-3)[Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc",
    expected = list(),
    line = 3631L
  ),
  list(
    pattern = ".-HexNAc$",
    glycan = "Neu5Ac(a2-3)Gal(b1-3)[Fuc(a1-4)]GlcNAc(b1-2)Man(a1-3)[Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc",
    expected = list("GlcNAc(b1-4)GlcNAc"),
    line = 3632L
  ),
  list(
    pattern = "!Neu5Ac-Gal-GlcNAc",
    glycan = "Neu5Ac(a2-3)Gal(b1-3)[Fuc(a1-4)]GlcNAc(b1-2)Man(a1-3)[Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc",
    expected = list("Gal(b1-4)GlcNAc"),
    line = 3633L
  ),
  list(
    pattern = "Hex-HexNAc-([Hex|Fuc]){1,2}-HexNAc(?=-HexNAc)",
    glycan = "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-3)[Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-4)GlcNAc",
    expected = list("Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc"),
    line = 3634L
  ),
  list(
    pattern = "Hex-HexNAc-([Hex|Fuc]){1,2}-HexNAc(?=-HexNAc)",
    glycan = "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-3)[Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-3)][Fuc(a1-6)]GlcNAc(b1-4)GlcNAc",
    expected = list("Man(b1-4)GlcNAc(b1-4)[Fuc(a1-3)][Fuc(a1-6)]GlcNAc"),
    line = 3635L
  ),
  list(
    pattern = "Hex-HexNAc-([Hex|Fuc]){1,2}-HexNAc(?!-HexNAc)",
    glycan = "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-3)[Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-3)][Fuc(a1-6)]GlcNAc",
    expected = list("Man(b1-4)GlcNAc(b1-4)[Fuc(a1-3)][Fuc(a1-6)]GlcNAc"),
    line = 3636L
  ),
  list(
    pattern = "(?<=Xyl-)Hex-HexNAc-([Hex|Fuc]){1,2}-HexNAc",
    glycan = "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-3)[Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-2)Man(a1-6)][Xyl(b1-2)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-3)]GlcNAc",
    expected = list("Man(b1-4)GlcNAc(b1-4)[Fuc(a1-3)]GlcNAc"),
    line = 3637L
  ),
  list(
    pattern = "(?<!Xyl-)Hex-HexNAc-([Hex|Fuc]){1,2}-HexNAc",
    glycan = "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-3)[Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-3)]GlcNAc",
    expected = list("Man(b1-4)GlcNAc(b1-4)[Fuc(a1-3)]GlcNAc"),
    line = 3638L
  ),
  list(
    pattern = "Hex%",
    glycan = "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc",
    expected = list("Gal", "Man", "Man"),
    line = 3640L
  ),
  list(
    pattern = "[HexNAc]{2}",
    glycan = "Gal(b1-4)GlcNAc(b1-4)GlcNAc",
    expected = list("GlcNAc(b1-4)GlcNAc"),
    line = 3642L
  ),
  list(
    pattern = "[HexNAc]{3}",
    glycan = "Gal(b1-4)GlcNAc(b1-4)GlcNAc",
    expected = list(),
    line = 3643L
  ),
  list(
    pattern = "[HexNAc]{2,}",
    glycan = "Gal(b1-4)GlcNAc(b1-4)GlcNAc(b1-4)GlcNAc",
    expected = list(
      "GlcNAc(b1-4)GlcNAc(b1-4)GlcNAc"
    ),
    line = 3645L
  ),
  list(
    pattern = "[Hex]{,2}",
    glycan = "Gal(b1-4)GlcNAc(b1-4)GlcNAc",
    expected = list("Gal"),
    line = 3646L
  ),
  list(
    pattern = "Hex-[HexNAc]{1,2}?",
    glycan = "Gal(b1-4)GlcNAc(b1-4)GlcNAc",
    expected = list(
      "Gal(b1-4)GlcNAc"
    ),
    line = 3648L
  ),
  list(
    pattern = "Hex-[HexNAc]*?",
    glycan = "Gal(b1-4)GlcNAc(b1-4)GlcNAc",
    expected = list(
      "Gal"
    ),
    line = 3649L
  ),
  list(
    pattern = "Hex-[HexNAc]+?",
    glycan = "Gal(b1-4)GlcNAc(b1-4)GlcNAc",
    expected = list(
      "Gal(b1-4)GlcNAc"
    ),
    line = 3650L
  ),
  list(
    pattern = "Mana6-Man",
    glycan = "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc",
    expected = list(
      "Man(a1-6)Man"
    ),
    line = 3652L
  ),
  list(
    pattern = "Mana3-Man",
    glycan = "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc",
    expected = list(
      "Man(a1-3)Man"
    ),
    line = 3653L
  ),
  list(
    pattern = "Galb4-GlcNAc",
    glycan = "Gal(b1-4)GlcNAc(b1-2)Man",
    expected = list(
      "Gal(b1-4)GlcNAc"
    ),
    line = 3654L
  ),
  list(
    pattern = "Neu5Aca3-Gal",
    glycan = "Neu5Ac(a2-3)Gal(b1-4)GlcNAc",
    expected = list(
      "Neu5Ac(a2-3)Gal"
    ),
    line = 3656L
  ),
  list(
    pattern = "Neu5Aca6-Gal",
    glycan = "Neu5Ac(a2-3)Gal(b1-4)GlcNAc",
    expected = list(),
    line = 3657L
  ),
  list(
    pattern = "[Gal|Man|Fuc]-HexNAc",
    glycan = "Gal(b1-4)GlcNAc",
    expected = list("Gal(b1-4)GlcNAc"),
    line = 3659L
  ),
  list(
    pattern = "[Gal|Man|Fuc]-HexNAc",
    glycan = "Fuc(a1-3)GlcNAc",
    expected = list("Fuc(a1-3)GlcNAc"),
    line = 3660L
  ),
  list(
    pattern = ".-.",
    glycan = "Gal(b1-4)GlcNAc",
    expected = list("Gal(b1-4)GlcNAc"),
    line = 3662L
  ),
  list(
    pattern = ".-.-.",
    glycan = "Gal(b1-4)GlcNAc",
    expected = list(),
    line = 3663L
  ),
  list(
    pattern = "^Hex-HexNAc$",
    glycan = "Gal(b1-4)GlcNAc",
    expected = TRUE,
    line = 3665L
  ),
  list(
    pattern = "^Hex-HexNAc$",
    glycan = "Gal(b1-4)GlcNAc(b1-4)GlcNAc",
    expected = FALSE,
    line = 3666L
  ),
  list(
    pattern = "Fuc-([^Gal])+-GlcNAc",
    glycan = "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-3)[Fuc(a1-3)[Fuc(a1-2)Gal(b1-4)]GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc",
    expected = list(),
    line = 3672L
  ),
  list(
    pattern = "Man-HexNAc",
    glycan = "Gal(b1-4)GlcNAc(b1-6)[Gal(b1-3)]GalNAc",
    expected = list(),
    line = 3673L
  ),
  list(
    pattern = "Hex-HexNAc-([Hex|Fuc]){1,2}-HexNAc(?=-HexNAc)",
    glycan = "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-3)[Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc",
    expected = list(),
    line = 3674L
  ),
  list(
    pattern = "Hex-HexNAc-([Hex|Fuc]){1,2}-HexNAc(?!-HexNAc)",
    glycan = "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-3)[Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-3)][Fuc(a1-6)]GlcNAc(b1-4)GlcNAc",
    expected = list(),
    line = 3675L
  ),
  list(
    pattern = "(?<=Xyl-)Hex-HexNAc-([Hex|Fuc]){1,2}-HexNAc",
    glycan = "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-3)[Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-3)]GlcNAc",
    expected = list(),
    line = 3676L
  ),
  list(
    pattern = "(?<!Xyl-)Hex-HexNAc-([Hex|Fuc]){1,2}-HexNAc",
    glycan = "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-3)[Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-2)Man(a1-6)][Xyl(b1-2)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-3)]GlcNAc",
    expected = list(),
    line = 3677L
  ),
  list(
    pattern = "Fuca3-([Galb4]){1}-GlcNAcb?",
    glycan = "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc",
    expected = list("Fuc(a1-3)[Gal(b1-4)]GlcNAc"),
    line = 3749L
  ),
  list(
    pattern = "Xylb2-([Mana3]){1}-([Mana6]){1}-Manb4-GlcNAcb4-GlcNAc",
    glycan = "Xyl(b1-2)[Man(a1-3)][Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc",
    expected = TRUE,
    line = 3762L
  ),
  list(
    pattern = "Galb3-([Sia(a2-6)]){1}-GalNAc",
    glycan = "Gal(b1-3)[Neu5Ac(a2-6)]GalNAc",
    expected = list(
      "Gal(b1-3)[Neu5Ac(a2-6)]GalNAc"
    ),
    line = 3765L
  ),
  list(
    pattern = "Fuc-([^Gal])+-GlcNAc",
    glycan = "Fuc(a1-3)[GlcNAc(b1-4)]GlcNAc",
    expected = list("Fuc(a1-3)[GlcNAc(b1-4)]GlcNAc"),
    line = 3786L
  ),
  list(
    pattern = "Fuc-([^Gal])+-GlcNAc",
    glycan = "Fuc(a1-3)[Gal(b1-4)]GlcNAc",
    expected = list(),
    line = 3787L
  ),
  list(
    pattern = "Fuc-([!Gal])+-GlcNAc",
    glycan = "Fuc(a1-3)[GlcNAc(b1-4)]GlcNAc",
    expected = list(
      "Fuc(a1-3)[GlcNAc(b1-4)]GlcNAc"
    ),
    line = 3788L
  ),
  list(
    pattern = "[HexNAc]{2,}",
    glycan = "Gal(b1-4)GlcNAc(b1-4)GlcNAc(b1-4)GlcNAc",
    expected = list("GlcNAc(b1-4)GlcNAc(b1-4)GlcNAc"),
    line = 3800L
  ),
  list(
    pattern = "Hex-HexNAc",
    glycan = "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc",
    expected = list("Gal(b1-4)GlcNAc", "Gal(b1-4)GlcNAc", "Man(b1-4)GlcNAc"),
    line = 3801L
  ),
  list(
    pattern = "Hex%",
    glycan = "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc",
    expected = list("Gal", "Man", "Man"),
    line = 3804L
  ),
  list(
    pattern = "Hex-[HexNAc]*?",
    glycan = "Gal(b1-4)GlcNAc(b1-4)GlcNAc",
    expected = list("Gal"),
    line = 3806L
  ),
  list(
    pattern = "Hex-[HexNAc]+?",
    glycan = "Gal(b1-4)GlcNAc(b1-4)GlcNAc",
    expected = list(
      "Gal(b1-4)GlcNAc"
    ),
    line = 3807L
  )
)
