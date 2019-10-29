# Code to diagnose the issue with joining bcicensusdat with nadja's growth+light dataframe
# Run at line 118 of workflow_june_newFGs.r

# Find mismatched tags
# Use 1995 data

cens1995 <- bcicensusdat[[3]] %>% select(tag, dbh)

cens_tag_integer <- as.integer(cens1995$tag)

cens_tag_1990_integer <- as.integer(bcicensusdat[[2]]$tag)

# Use set ops to see where the mismatches are


length(setdiff(cens_tag_integer, growth9095$tag))

length(intersect(cens_tag_integer, growth9095$tag))



length(intersect(cens_tag_integer, growth8590$tag))

newrecruits1995 <- setdiff(cens_tag_integer, cens_tag_1990_integer)

table(newrecruits1995 %in% growth9095$tag)
