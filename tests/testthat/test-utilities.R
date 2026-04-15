# =============================================================================
# Tests: utility functies
# validate_parentheses, clean_constraints, expand_parentheses,
# process_abs_and_expand, process_constraint_syntax,
# expand_compound_constraints, list_to_df_rows, format_numeric,
# detect_range_restrictions, remove_redundant_constraints,
# check_sample_nobs, check_N_with_sample_nobs, tukeyChi
# =============================================================================


# =============================================================================
# validate_parentheses
# =============================================================================

test_that("validate_parentheses: gebalanceerde haakjes", {
  expect_true(restriktor:::validate_parentheses("(a, b)"))
  expect_true(restriktor:::validate_parentheses("((a))"))
  expect_true(restriktor:::validate_parentheses("a + b"))
  expect_true(restriktor:::validate_parentheses(""))
})

test_that("validate_parentheses: ongebalanceerde haakjes", {
  expect_false(restriktor:::validate_parentheses("(a"))
  expect_false(restriktor:::validate_parentheses("a)"))
  expect_false(restriktor:::validate_parentheses(")("))
  expect_false(restriktor:::validate_parentheses("((a)"))
})

test_that("validate_parentheses: geneste haakjes", {
  expect_true(restriktor:::validate_parentheses("((a, b), (c, d))"))
  expect_false(restriktor:::validate_parentheses("((a, b), (c, d)"))
})


# =============================================================================
# clean_constraints
# =============================================================================

test_that("clean_constraints: verwijdert comments", {
  result <- restriktor:::clean_constraints("x1>x2 # dit is een opmerking\nx3>x4")
  expect_false(grepl("#", result))
  expect_false(grepl("opmerking", result))
})

test_that("clean_constraints: normaliseert delimiters", {
  result <- restriktor:::clean_constraints("x1>x2; x3>x4")
  expect_true(grepl("\n", result))
  expect_false(grepl(";", result))

  result2 <- restriktor:::clean_constraints("x1>x2 & x3>x4")
  expect_false(grepl("&", result2))
})

test_that("clean_constraints: verwijdert spaties", {
  result <- restriktor:::clean_constraints("x1 > x2")
  expect_false(grepl(" ", result))
  expect_true(grepl("x1>x2", result))
})

test_that("clean_constraints: komma buiten haakjes wordt newline", {
  result <- restriktor:::clean_constraints("x1>x2, x3>x4")
  expect_true(grepl("\n", result))
})

test_that("clean_constraints: komma binnen haakjes blijft behouden", {
  result <- restriktor:::clean_constraints("(x1, x2) > x3")
  expect_true(grepl(",", result))
})

test_that("clean_constraints: comment op laatste regel wordt verwijderd", {
  # Regressie: (?=\n) matcht niet op einde-van-string.
  # Fix: (?=\n|$)
  result <- restriktor:::clean_constraints("x1>x2 # comment op laatste regel")
  expect_false(grepl("#", result))
  expect_false(grepl("comment", result))
  expect_true(grepl("x1>x2", result))
})


# =============================================================================
# expand_parentheses
# =============================================================================

test_that("expand_parentheses: breidt eenvoudige haakjes uit", {
  result <- restriktor:::expand_parentheses("(x1,x2)>x3")
  expect_equal(length(result), 2)
  expect_true("x1>x3" %in% result)
  expect_true("x2>x3" %in% result)
})

test_that("expand_parentheses: zonder haakjes geeft origineel terug", {
  result <- restriktor:::expand_parentheses("x1>x2")
  expect_equal(result, "x1>x2")
})

test_that("expand_parentheses: fout bij ongebalanceerde haakjes", {
  expect_error(
    restriktor:::expand_parentheses("(x1,x2>x3"),
    "parenthes"
  )
})

test_that("expand_parentheses: geneste haakjes", {
  result <- restriktor:::expand_parentheses("(x1,x2)>(x3,x4)")
  expect_equal(length(result), 4)
})


# --- Regressietest: rekenkundige haakjes (fix expand_parentheses) -----------
# Bug: (A - B) / 1.516 > 0.2 werd foutief uitgebreid naar A - B / 1.516 > 0.2
# waardoor deling alleen op B werd toegepast, niet op A.
# Fix: haakjes zonder komma = rekenkundige groepering -> niet uitbreiden.

test_that("expand_parentheses: rekenkundige haakjes blijven intact (regressie)", {
  # Kerncase: deling na haakjes moet NIET de haakjes strippen
  result <- restriktor:::expand_parentheses("(GroupPassive-GroupActive)/1.516>0.2")
  expect_equal(result, "(GroupPassive-GroupActive)/1.516>0.2")
  expect_equal(length(result), 1)
})

test_that("expand_parentheses: vermenigvuldiging na haakjes blijft intact", {
  result <- restriktor:::expand_parentheses("2*(x1-x2)>0")
  expect_equal(result, "2*(x1-x2)>0")
})

test_that("expand_parentheses: geneste rekenkundige haakjes", {
  result <- restriktor:::expand_parentheses("((x1-x2)+(x3-x4))/2>0")
  expect_equal(result, "((x1-x2)+(x3-x4))/2>0")
})

test_that("expand_parentheses: redundante haakjes zonder operator", {
  # (x1-x2)>0 is ook rekenkundig, geen komma -> intact
  result <- restriktor:::expand_parentheses("(x1-x2)>0")
  expect_equal(result, "(x1-x2)>0")
})

test_that("expand_parentheses: enkel geneste rekenkundige haakjes", {
  # Dubbel-geneste haakjes ((a,b),(c,d)) worden niet uitgebreid
  # omdat de functie het eerste en tweede haakje pakt (niet matching paar).
  # Dit is een bekende beperking. Gebruik in dit geval apart:
  #   (x1-x2) > 0
  #   (x3-x4) > 0
  result <- restriktor:::expand_parentheses("((x1-x2),(x3-x4))>0")
  # Retourneert ongewijzigd vanwege geneste haakjes beperking
  expect_equal(result, "((x1-x2),(x3-x4))>0")
})

test_that("expand_parentheses: komma-expansie werkt nog steeds correct", {
  # Controleer dat bain-stijl niet gebroken is door de fix
  result <- restriktor:::expand_parentheses("(x1,x2,x3)>0")
  expect_equal(length(result), 3)
  expect_true("x1>0" %in% result)
  expect_true("x2>0" %in% result)
  expect_true("x3>0" %in% result)
})


# =============================================================================
# process_abs_and_expand
# =============================================================================

test_that("process_abs_and_expand: abs() wordt beschermd tegen expansie", {
  result <- restriktor:::process_abs_and_expand("abs(x1-x2)>0")
  expect_true(any(grepl("abs\\(x1-x2\\)", result)))
})

test_that("process_abs_and_expand: zonder abs() werkt normaal", {
  result <- restriktor:::process_abs_and_expand("(x1,x2)>x3")
  expect_equal(length(result), 2)
})

test_that("process_abs_and_expand: abs() met haakjes-expansie", {
  result <- restriktor:::process_abs_and_expand("abs(x1-x2)>(x3,x4)")
  expect_equal(length(result), 2)
  expect_true(all(grepl("abs\\(x1-x2\\)", result)))
})

test_that("process_abs_and_expand: abs() met geneste haakjes (regressie)", {
  # Regressie: regex [^)]+ stopt bij eerste ')' in abs(x1-(x2+x3))
  # Fix: recursieve regex voor gebalanceerde haakjes
  result <- restriktor:::process_abs_and_expand("abs(x1-(x2+x3))>0")
  expect_true(any(grepl("abs\\(x1-\\(x2\\+x3\\)\\)", result)))
  expect_equal(length(result), 1)
})


# =============================================================================
# process_constraint_syntax
# =============================================================================

test_that("process_constraint_syntax: normaliseert == naar =", {
  result <- restriktor:::process_constraint_syntax("x1==x2")
  expect_true(any(grepl("x1=x2", result)))
  expect_false(any(grepl("==", result)))
})

test_that("process_constraint_syntax: normaliseert <= naar <", {
  result <- restriktor:::process_constraint_syntax("x1<=x2")
  expect_true(any(grepl("x1<x2", result)))
  expect_false(any(grepl("<=", result)))
})

test_that("process_constraint_syntax: normaliseert >= naar >", {
  result <- restriktor:::process_constraint_syntax("x1>=x2")
  expect_true(any(grepl("x1>x2", result)))
  expect_false(any(grepl(">=", result)))
})

test_that("process_constraint_syntax: verwijdert lege regels", {
  result <- restriktor:::process_constraint_syntax("x1>x2\n\nx3>x4")
  expect_false(any(result == ""))
  expect_equal(length(result), 2)
})

test_that("process_constraint_syntax: split op newline", {
  result <- restriktor:::process_constraint_syntax("x1>x2\nx3=x4\nx5<x6")
  expect_equal(length(result), 3)
})


# =============================================================================
# expand_compound_constraints
# =============================================================================

test_that("expand_compound_constraints: ketenvergelijking x1>x2>x3", {
  result <- restriktor:::expand_compound_constraints(list("x1>x2>x3"))
  expect_true(is.list(result))
  expanded <- result[[1]]
  expect_equal(length(expanded), 2)
  expect_true(any(grepl("x1>x2", expanded)))
  expect_true(any(grepl("x2>x3", expanded)))
})

test_that("expand_compound_constraints: enkele vergelijking blijft intact", {
  result <- restriktor:::expand_compound_constraints(list("x1>x2"))
  expect_equal(result[[1]], "x1>x2")
})

test_that("expand_compound_constraints: gelijkheidsketen x1=x2=x3", {
  result <- restriktor:::expand_compound_constraints(list("x1=x2=x3"))
  expanded <- result[[1]]
  expect_equal(length(expanded), 2)
})

test_that("expand_compound_constraints: meerdere hypothesen tegelijk", {
  result <- restriktor:::expand_compound_constraints(list("x1>x2>x3", "a=b"))
  expect_equal(length(result), 2)
  expect_equal(length(result[[1]]), 2)
  expect_equal(result[[2]], "a=b")
})


# =============================================================================
# list_to_df_rows
# =============================================================================

test_that("list_to_df_rows: basisconversie", {
  lst <- list(a = c(x = 1, y = 2), b = c(x = 3, y = 4))
  result <- restriktor:::list_to_df_rows(lst)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
})

test_that("list_to_df_rows: ongelijke namen worden samengevoegd", {
  lst <- list(a = c(x = 1), b = c(x = 2, y = 3))
  result <- restriktor:::list_to_df_rows(lst)
  expect_equal(ncol(result), 2)
  expect_true(is.na(result[1, "y"]))
})

test_that("list_to_df_rows: fout bij lege lijst", {
  expect_error(restriktor:::list_to_df_rows(list()))
})


# =============================================================================
# format_numeric
# =============================================================================

test_that("format_numeric: standaard formatting", {
  result <- restriktor:::format_numeric(3.14159, digits = 3)
  expect_true(grepl("3.142", result))
})

test_that("format_numeric: nul wordt correct geformatted", {
  result <- restriktor:::format_numeric(0, digits = 3)
  expect_true(grepl("0", result))
})

test_that("format_numeric: zeer kleine waarden in scientific notatie", {
  result <- restriktor:::format_numeric(0.00001, digits = 3)
  expect_true(grepl("[eE]", result))
})

test_that("format_numeric: grote waarden in scientific notatie", {
  result <- restriktor:::format_numeric(50000, digits = 3)
  expect_true(grepl("[eE]", result))
})

test_that("format_numeric: near-zero wordt 0", {
  result <- restriktor:::format_numeric(1e-10, digits = 3)
  expect_true(grepl("0", result))
})


# =============================================================================
# detect_range_restrictions
# =============================================================================

test_that("detect_range_restrictions: detecteert range (rij i = -rij j)", {
  # x1 >= 0 en x1 <= 1 => rij1 = c(1,0), rij2 = c(-1,0)
  Amat <- rbind(c(1, 0), c(-1, 0))
  result <- restriktor:::detect_range_restrictions(Amat)
  expect_true(nrow(result) > 0)
  expect_equal(result[1, ], c(1, 2))
})

test_that("detect_range_restrictions: geen range geeft lege matrix", {
  Amat <- rbind(c(1, 0), c(0, 1))
  result <- restriktor:::detect_range_restrictions(Amat)
  expect_equal(nrow(result), 0)
})

test_that("detect_range_restrictions: meerdere ranges", {
  # x1 in [0,1] en x2 in [0,1]
  Amat <- rbind(c(1, 0), c(-1, 0), c(0, 1), c(0, -1))
  result <- restriktor:::detect_range_restrictions(Amat)
  expect_equal(nrow(result), 2)
})


# =============================================================================
# remove_redundant_constraints
# =============================================================================

test_that("remove_redundant_constraints: verwijdert redundante rijen", {
  # x1 >= 1 en x1 >= 2 => alleen x1 >= 2 is nodig (strictst)
  constraints <- rbind(c(1, 0), c(1, 0))
  rhs <- c(1, 2)
  result <- restriktor:::remove_redundant_constraints(constraints, rhs)
  expect_equal(nrow(result$constraints), 1)
  expect_equal(result$rhs, 2)
})

test_that("remove_redundant_constraints: behoudt unieke restricties", {
  constraints <- rbind(c(1, 0), c(0, 1))
  rhs <- c(0, 0)
  result <- restriktor:::remove_redundant_constraints(constraints, rhs)
  expect_equal(nrow(result$constraints), 2)
})

test_that("remove_redundant_constraints: drie redundante restricties", {
  constraints <- rbind(c(1, 0), c(1, 0), c(1, 0))
  rhs <- c(1, 2, 3)
  result <- restriktor:::remove_redundant_constraints(constraints, rhs)
  expect_equal(nrow(result$constraints), 1)
  # Hoogste rhs (meest restrictieve) wordt behouden
  expect_equal(result$rhs, 3)
})


# =============================================================================
# check_sample_nobs
# =============================================================================

test_that("check_sample_nobs: scalar wordt doorgegeven", {
  result <- restriktor:::check_sample_nobs(100)
  expect_equal(result, 100)
})

test_that("check_sample_nobs: vector wordt gesommeerd", {
  expect_message(
    result <- restriktor:::check_sample_nobs(c(30, 30, 30)),
    "group sizes"
  )
  expect_equal(result, 90)
})

test_that("check_N_with_sample_nobs: gelijke waarden geven N terug", {
  result <- restriktor:::check_N_with_sample_nobs(100, 100)
  expect_equal(result, 100)
})

test_that("check_N_with_sample_nobs: ongelijke waarden geven N terug met message", {
  expect_message(
    result <- restriktor:::check_N_with_sample_nobs(100, 50),
    "differs"
  )
  expect_equal(result, 100)
})

test_that("check_N_with_sample_nobs: NULL sample_nobs geeft N terug", {
  result <- restriktor:::check_N_with_sample_nobs(100, NULL)
  expect_equal(result, 100)
})


# =============================================================================
# tukeyChi
# =============================================================================

test_that("tukeyChi: rho functie (deriv=0) bereik [0,1]", {
  x <- seq(-10, 10, by = 0.5)
  rho <- restriktor:::tukeyChi(x, deriv = 0)
  expect_true(all(rho >= 0 - 1e-10))
  expect_true(all(rho <= 1 + 1e-10))
})

test_that("tukeyChi: rho=1 voor |x| > c", {
  rho <- restriktor:::tukeyChi(10, c = 4.685, deriv = 0)
  expect_equal(rho, 1)
})

test_that("tukeyChi: rho=0 bij x=0", {
  rho <- restriktor:::tukeyChi(0, deriv = 0)
  expect_equal(rho, 0)
})

test_that("tukeyChi: psi functie (deriv=1) is 0 buiten [-c,c]", {
  psi <- restriktor:::tukeyChi(10, c = 4.685, deriv = 1)
  expect_equal(psi, 0)
})

test_that("tukeyChi: fout bij ongeldige deriv", {
  expect_error(restriktor:::tukeyChi(1, deriv = 3), "deriv")
})
