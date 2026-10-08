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

