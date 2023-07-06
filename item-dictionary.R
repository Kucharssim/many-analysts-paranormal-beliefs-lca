item_dictionary <- list(
  c("par_1",  0, 2, 0, 0, 2, "Predictive dreams"),
  c("par_2",  0, 2, 1, 0, 2, "Contact with deceased"),
  c("par_3",  1, 2, 1, 0, 2, "Astrology"),
  c("par_4",  0, 2, 2, 0, 3, "Spirits and ghosts"),
  c("par_5",  2, 0, 0, 0, 1, "Healing by laying on of hands"),
  c("par_6",  1, 2, 1, 0, 2, "Palmistry"),
  c("par_7",  2, 1, 0, 0, 1, "Seeing into the future"),
  c("par_8",  2, 0, 0, 0, 1, "Life course description"),
  c("par_9",  2, 0, 0, 0, 1, "Mind-reading"),
  c("par_10", 2, 1, 0, 0, 1, "Seeing events happening"),
  c("par_11", 2, 0, 0, 0, 1, "Graphology"),
  c("par_12", 0, 0, 2, 0, 3, "Extraterrestrial beings"),
  c("par_13", 2, 0, 0, 0, 1, "Dowsing"),
  c("par_14", NA, NA, NA, NA, 2, "Telekinesis"),
  c("par_15", 2, 1, 1, 0, 1, "Causing events (”willpower”)"),
  c("par_16", 0, 1, 2, 0, 3, "Gnomes and elves"),
  c("par_17", 0, 2, 1, 0, 2, "Reincarnation"),
  c("par_18", 1, 0, 2, 0, 3, "Sorcery, black magic"),
  c("par_19", 0, 1, 2, 0, 3, "Haunted places"),
  c("par_20", 0, 0, 0, 2, 4, "Lucky numbers"),
  c("par_21", 0, 0, 0, 2, 4, "Unlucky numbers"),
  c("par_22", 0, 0, 1, 2, 4, "Walking under a ladder"),
  c("par_23", 0, 2, 1, 0, 2, "Recounting (previous life)"),
  c("par_24", 0, 0, 2, 0, 3, "Devil possesion")
)
item_dictionary <- do.call(rbind, item_dictionary) |> as.data.frame()
colnames(item_dictionary) <- c("item", "eha", "sr", "ub", "es", "factor", "label")
item_dictionary$item <- factor(item_dictionary$item, levels = sprintf("par_%s", 1:24))
item_dictionary$eha  <- as.numeric(item_dictionary$eha)
item_dictionary$sr   <- as.numeric(item_dictionary$sr)
item_dictionary$ub   <- as.numeric(item_dictionary$ub)
item_dictionary$es   <- as.numeric(item_dictionary$es)
item_dictionary$factor <- factor(
  item_dictionary$factor,
  levels = 1:4,
  labels = c("Extraordinary human abilities", "Supernatural reality", "Unearthly beings", "Everyday superstition")
  )
