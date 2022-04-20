common_to_binomial <- function(x){
  species_dictionary <- list(
    armadillo = "Dasypus novemcinctus",
    bobcat = "Lynx rufus",
    eastern_chipmunk = "Tamias striatus",
    california_ground_squirrel = "Otospermophilus beecheyi",
    cottontail_sp = c("Sylvilagus bachmani",
                      "Sylvilagus floridanus",
                      "Sylvilagus audubonii",
                      "Sylvilagus palustris"),
    coyote = "Canis latrans",
    gray_squirrel_sp = c("Sciurus carolinensis",
                         "Sciurus griseus"),
    gray_fox = "Urocyon cinereoargenteus",
    weasel_sp = c("Mustela nivalis",
                  "Mustela frenata"),
    muskrat = "Ondatra zibethicus",
    north_american_beaver = "Castor canadensis",
    north_american_river_otter = "Lontra canadensis",
    raccoon = "Procyon lotor",
    red_fox = "Vulpes vulpes",
    virginia_opossum = "Didelphis virginiana",
    white_tailed_deer = "Odocoileus virginianus",
    woodchuck = "Marmota monax",
    jackrabbit_sp = c("Lepus townsendii",
                      "Lepus californicus",
                      "Lepus alleni",
                      "Lepus americanus"
                      ),
    feral_hog = "Sus scrofa",
    fox_squirrel = "Sciurus niger",
    ringtail = "Bassariscus astutus",
    striped_skunk = "Mephitis mephitis",
    mule_deer = "Odocoileus hemionus",
    red_squirrel = "Tamiasciurus hudsonicus",
    flying_squirrel_sp = c("Glaucomys sabrinus",
                           "Glaucomys volans",
                           "Glaucomys oregonensis"),
    north_american_mink = "Neovison vison",
    black_bear = "Ursus americanus",
    elk = "Cervus canadensis",
    cougar = "Puma concolor",
    black_tailed_prairie_dog = "Cynomys ludovicianus",
    least_chipmunk = "Neotamias minimus",
    moose = "Alces alces",
    north_american_porcupine = "Erethizon dorsatum",
    richardson_ground_squirrel = "Urocitellus richardsonii",
    north_american_badger = "Taxidea taxus",
    thirteen_lined_ground_squirrel = "Ictidomys tridecemlineatus",
    black_tailed_deer = "Odocoileus hemionus",
    american_hog_nosed_skunk = "Conepatus leuconotus",
    harris_antelope_squirrel = "Ammospermophilus harrisii",
    hooded_skunk = "Mephitis macroura",
    javelina = "Pecari tajacu",
    kit_fox = "Vulpes macrotis",
    pronghorn = "Antilocapra americana",
    rock_squirrel = "Otospermophilus variegatus",
    round_tailed_ground_squirrel = "Xerospermophilus tereticaudus",
    western_spotted_skunk = "Spilogale gracilis",
    white_nosed_coatimundi = "Nasua narica",
    fisher = "Martes pennanti",
    douglas_squirrel = "Tamiasciurus douglasii",
    nutria = "Myocastor coypus"
  )
  standardized_names <- names(species_dictionary)
  
  x <- unique(x)
  
  species_location <- which(standardized_names %in% x)
  
  to_return <- utils::stack(species_dictionary[species_location])
  colnames(to_return) <- c("Binomial", "Species")
  return(to_return)
  
}
