# Unit Tests

Test case for the deduper

- `deduper_test_input.sam`: sorted input, 11 reads
- `deduper_test_output.sam`: expected output after deduplication, 6 reads

Every read comes from the provided `test.sam`. The header is copied from it, with only the sort tag changed and the fields needed for a test case changed. 

## Test cases

QNAMEs shortened to the last three fields. All on chromosome 2 unless noted.

| # | QNAME (end) | Changed | Strand | POS | CIGAR | 5' start | Result | Case tested |
|---|---|---|---|---|---|---|---|---|
| 1 | 25071:87858:TAGCAAGG | nothing | + | 76718924 | 71M | 76718924 | keep | first of its kind |
| 2 | 38067:55266:GAACAGGT | UMI | + | 76718924 | 71M | 76718924 | keep | same position as 1, different UMI so not duplicate |
| 3 | 9957:51931:TAGCAAGG | FLAG to 16 | - | 76718924 | 71M | 76718924+71-1 = 76718994 | keep | same position and UMI as read 1, different strand, can't be duplicate |
| 4 | 52657:6101:TAGCAAGG | FLAG to 16, CIGAR, MD | - | 76718924 | 69M2S | 76718924+69+2-1 = 76718994 | remove | right soft clip on minus strand, matches 3 when accounting for soft clip/5' start|
| 5 | 4191:1194:TAGCAAGG | POS, CIGAR, MD | + | 76718927 | 3S68M | 76718927-3 = 76718924 | remove | left soft clip means that 5' position is required not first aligned base. without 5' factored can't see that it matches read 1 with some lack of alignment/possible sequence error |
| 6 | 94095:71756:AACGCCAT | nothing | + | 76875967 | 15M470N56M | 76875967 | keep | No earlier read with this position, strand, and UMI |
| 7 | 69992:67325:AACGCCAT | nothing | + | 76875967 | 15M470N56M | 76875967 | remove | real duplicate |
| 8 | 347:19169:AACGCCAT | FLAG to 16 | - | 76875967 | 15M470N56M | 76875967+15+470+56-1 = 76876507 | keep | Same position and UMI as read 6, different strand, can’t be duplicate |
| 9 | 6610:1112:AACGCCAT | FLAG to 16, POS, CIGAR, MD | - | 76875969 | 2S13M470N56M | 76875969+13+470+56-1 = 76876507 | remove | left soft clip doesn’t count on the minus strand, its right end matches read 8’s along with other info so duplicate |
| 10 | 5571:1221:GTGATGTC | nothing | + | 76901561 | 71M | 76901561 | keep | first of its kind |
| 11 | 20566:1080:TATAGNGC | nothing | + | 119767599 | 71M | n/a | discard | unknown UMI (N in UMI) |
| 12 | 18996:1145:TTCGCCTA | nothing | + | 130171653 | 40M1I30M | 130171653 | keep | Has an insertion, which doesn’t change a plus-strand 5' start, tests can handle insertions. No earlier read matches so is kept |
| 13 | 86637:67646:GTGATGTC | RNAME to 3 | + | 76901561 | 71M | 76901561 | keep | Same position, strand, and UMI as read 10, different chromosome, can’t be duplicate |

## End Result

- 13 reads in
- 8 kept (1, 2, 3, 6, 8, 10, 12, 13)
- 4 duplicates removed (4, 5, 7, 9)
- 1 discarded for unknown UMI (11)

