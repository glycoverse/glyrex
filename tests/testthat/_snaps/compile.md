# malformed patterns fail at compilation

    Code
      rex_compile("[Gal")
    Condition
      Error:
      ! Invalid pattern at character 5: Expected ']'.

---

    Code
      rex_compile("[Gal]{3,2}")
    Condition
      Error:
      ! Invalid pattern at character 11: Invalid repetition bounds.

---

    Code
      rex_compile("[Gal|]")
    Condition
      Error:
      ! Invalid pattern at character 6: Empty pattern or alternative.

---

    Code
      rex_compile("Gal*")
    Condition
      Error:
      ! Invalid pattern at character 5: Attach repetition to a bracketed group.

---

    Code
      rex_compile("(?<x>Gal)-(?<x>Man)")
    Condition
      Error:
      ! Invalid pattern at character 15: Invalid, duplicate, or reserved capture name.

---

    Code
      rex_compile("[GlcNAc]{,}")
    Condition
      Error:
      ! Invalid pattern at character 12: Invalid repetition range.

---

    Code
      rex_compile("(?=)")
    Condition
      Error:
      ! Invalid pattern at character 4: Empty pattern or alternative.

---

    Code
      rex_compile(NA_character_)
    Condition
      Error:
      ! `pattern` must be one non-empty, non-missing string.

# the full-match column name cannot be shadowed by a capture

    Code
      rex_compile("(?<match>Gal)")
    Condition
      Error:
      ! Invalid pattern at character 9: Invalid, duplicate, or reserved capture name.

# bracketed exclamation negation is rejected

    Code
      rex_compile("[!Fuc]")
    Condition
      Error:
      ! Invalid pattern at character 2: '[!...]' is not supported; use '[^...]' for residue negation or '(!...-)' for branch absence.

---

    Code
      rex_compile("Gal-([!Fuca3])-GlcNAc")
    Condition
      Error:
      ! Invalid pattern at character 7: '[!...]' is not supported; use '[^...]' for residue negation or '(!...-)' for branch absence.

---

    Code
      rex_compile("[Gal|!Man]")
    Condition
      Error:
      ! Invalid pattern at character 6: '[!...]' is not supported; use '[^...]' for residue negation or '(!...-)' for branch absence.

# only explicit branch absence groups are accepted

    Code
      rex_compile("Gal-!Fuca3-GlcNAc")
    Condition
      Error:
      ! Invalid pattern at character 5: Bare !residue is not supported; use (!branch-) for branch absence or [^residue] for residue negation.

---

    Code
      rex_compile("!Fuc")
    Condition
      Error:
      ! Invalid pattern at character 1: Bare !residue is not supported; use (!branch-) for branch absence or [^residue] for residue negation.

---

    Code
      rex_compile("Gal-(!Fuc)-GlcNAc")
    Condition
      Error:
      ! Invalid pattern at character 10: Branch absence requires a trailing attachment dash: use (!branch-)residue.

---

    Code
      rex_compile("Gal-(!Fuca3)-GlcNAc")
    Condition
      Error:
      ! Invalid pattern at character 12: Branch absence requires a trailing attachment dash: use (!branch-)residue.

# patterns reject IUPAC-style linkages

    Code
      rex_compile("Gal(b1-4)GlcNAc")
    Condition
      Error:
      ! Invalid pattern at character 4: IUPAC-style linkages are not supported in patterns; use glyco-regex syntax such as 'Galb4-GlcNAc' or 'GlcNAca$'.

---

    Code
      rex_compile("Gal(b1-4)-GlcNAc")
    Condition
      Error:
      ! Invalid pattern at character 4: IUPAC-style linkages are not supported in patterns; use glyco-regex syntax such as 'Galb4-GlcNAc' or 'GlcNAca$'.

---

    Code
      rex_compile("GlcNAc(a1-")
    Condition
      Error:
      ! Invalid pattern at character 7: IUPAC-style linkages are not supported in patterns; use glyco-regex syntax such as 'Galb4-GlcNAc' or 'GlcNAca$'.

---

    Code
      rex_compile("[GlcNAc(b1-)]")
    Condition
      Error:
      ! Invalid pattern at character 8: IUPAC-style linkages are not supported in patterns; use glyco-regex syntax such as 'Galb4-GlcNAc' or 'GlcNAca$'.

---

    Code
      rex_compile("Gal-(!Fuc(a1-3)-)GlcNAc")
    Condition
      Error:
      ! Invalid pattern at character 10: IUPAC-style linkages are not supported in patterns; use glyco-regex syntax such as 'Galb4-GlcNAc' or 'GlcNAca$'.

---

    Code
      rex_compile("[Gal(?1-3/4)]-GlcNAc")
    Condition
      Error:
      ! Invalid pattern at character 5: IUPAC-style linkages are not supported in patterns; use glyco-regex syntax such as 'Galb4-GlcNAc' or 'GlcNAca$'.

---

    Code
      rex_detect("Gal(b1-4)GlcNAc", "Gal(b1-4)GlcNAc")
    Condition
      Error:
      ! Invalid pattern at character 4: IUPAC-style linkages are not supported in patterns; use glyco-regex syntax such as 'Galb4-GlcNAc' or 'GlcNAca$'.

# branches require the attachment dash inside parentheses

    Code
      rex_compile("Galb4-(Fuca3)-GlcNAc")
    Condition
      Error:
      ! Invalid pattern at character 13: Branches require a trailing attachment dash: use (branch-)residue.

---

    Code
      rex_compile("Galb4-(Fuca3)GlcNAc")
    Condition
      Error:
      ! Invalid pattern at character 13: Branches require a trailing attachment dash: use (branch-)residue.

---

    Code
      rex_compile("Galb4-(Fuca3-)-GlcNAc")
    Condition
      Error:
      ! Invalid pattern at character 15: Put the branch attachment dash inside parentheses: use (branch-)residue.

---

    Code
      rex_compile("Galb4-(Fuca3-){1}-GlcNAc")
    Condition
      Error:
      ! Invalid pattern at character 18: Put the branch attachment dash inside parentheses: use (branch-)residue.

