
### ensemble restart files
res_file1_info = {
    'filepath': '/data01/shzhu/tagch4_lf/Restart/',
    'filename': 'Restart.ST000.EN0001-EN0047.YYYYMMDD_0000z.nc4',
    'var_prefix': ['SpeciesRst_CH4_TAG_AN','SpeciesRst_CH4_TAG_NA'],
    'old_item': [1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,\
                 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47],
    'new_item': [1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,\
                 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47]
}

res_file2_info = {
    'filepath': '/data01/shzhu/tagch4_lf/Restart/',
    'filename': 'Restart.ST000.EN0048-EN0093.YYYYMMDD_0000z.nc4',
    'var_prefix': ['SpeciesRst_CH4_TAG_AN','SpeciesRst_CH4_TAG_NA'],
    'old_item': [1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,\
                 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46],
    'new_item': [48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71,\
                 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93]
}

res_info_list = [res_file1_info, res_file2_info]
# res_file1_info = {
#     'filepath': '/data/shzhu/tagch4/Restart/',
#     'filename': 'Restart.ST000.EN0001-EN0037.YYYYMMDD_0000z.nc4',
#     'var_prefix': ['SpeciesRst_CH4_TAG_AN'],
#     'old_item': [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31],
#     'new_item': [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31]
# }

# res_file3_info = {
#     'filepath': '/data/shzhu/tagch4_2x25_AN/Restart/',
#     'filename': 'Restart.ST000.EN0001-EN0037.YYYYMMDD_0000z.nc4',
#     'var_prefix': ['SpeciesRst_CH4_TAG_NA'],
#     'old_item': [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31],
#     'new_item': [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31]
# }
# res_file2_info = {
#     'filepath': '/data/shzhu/tagch4/Restart/',
#     'filename': 'Restart.ST000.EN0038-EN0082.YYYYMMDD_0000z.nc4',
#     'var_prefix': ['SpeciesRst_CH4_TAG_AN'],
#     'old_item': [ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,\
#                   23, 24, 29, 30, 31, 32, 33, 34, 35, 37, 38, 39, 40, 41, 42, 43, 44, 45],
#     'new_item': [48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69,\
#                  70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81,82, 83, 84, 85, 86, 87]
# }

# res_file4_info = {
#     'filepath': '/data/shzhu/tagch4_2x25_AN/Restart/',
#     'filename': 'Restart.ST000.EN0038-EN0082.YYYYMMDD_0000z.nc4',
#     'var_prefix': ['SpeciesRst_CH4_TAG_NA'],
#     'old_item':  [ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,\
#                   23, 24, 29, 30, 31, 32, 33, 34, 35, 37, 38, 39, 40, 41, 42, 43, 44, 45],
#     'new_item': [48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69,\
#                  70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87]
# }


# res_file5_info = {
#     'filepath': '/data/shzhu/tagch4_new/Restart/',
#     'filename': 'Restart.ST000.EN0001-EN0022.YYYYMMDD_0000z.nc4',
#     'var_prefix': ['SpeciesRst_CH4_TAG_AN','SpeciesRst_CH4_TAG_NA'],
#     'old_item': [1, 2, 3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22],
#     'new_item': [32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 88, 89, 90, 91, 92, 93]
# }

# res_info_list = [res_file1_info, res_file2_info, res_file3_info, res_file4_info, res_file5_info]