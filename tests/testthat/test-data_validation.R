
# Start tests ----
test_that("check_geneid_match() flags mismatches between gene sets", {
    
    set1 <- c("gene1", "gene2A", "gene3", "gene4A")
    set2 <- c("gene1", "gene2", "gene3", "gene4")
    
    expect_error(check_geneid_match(set1, set2))
})
