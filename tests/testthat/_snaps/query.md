# unsupported inputs and incompatible recycling fail clearly

    Code
      rex_detect(c("Gal", "Man"), c("Gal", "Man", "Fuc"))
    Condition
      Error:
      ! Inputs must have equal lengths or length one.

---

    Code
      rex_subset("Gal", c("Gal", "Man"))
    Condition
      Error:
      ! Subsetting requires one pattern or one per glycan.

---

    Code
      rex_detect(1, "Gal")
    Condition
      Error:
      ! `string` must be character, glyrepr_structure, or a glycan igraph.

---

    Code
      rex_detect("Gal", "Gal", negate = NA)
    Condition
      Error:
      ! `negate` must be TRUE or FALSE.

---

    Code
      rex_detect("{Fuc(a1-3)}Gal(b1-4)GlcNAc", "Gal")
    Condition
      Error:
      ! Queries require a connected rooted tree without floating metadata.

# floating substituents are rejected for every representation

    Code
      rex_detect(x, "Gal")
    Condition
      Error:
      ! Queries require a connected rooted tree without floating metadata.

---

    Code
      rex_count(as.list(x)[[1]], "Gal")
    Condition
      Error:
      ! Queries require a connected rooted tree without floating metadata.

---

    Code
      rex_extract(as.character(x), "Gal")
    Condition
      Error:
      ! Queries require a connected rooted tree without floating metadata.

