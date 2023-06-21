
include {
    CygnusConsensus
 } from "./modules/cygnus"


include {
    RotateBySequence
} from "../parse_convert/modules/rotators"

workflow Cygnus {
    take:
        reads
        adapters
    main:
        cygnus_consensus = CygnusConsensus(reads)
        RotateBySequence(cygnus_consensus)

    emit:
        cygnus_consensus

}