import sys

from hetdex_api.survey import Survey

HDRVersion = "hdr5"

if __name__ == '__main__':

    min_month = int(sys.argv[1])
    max_month = int(sys.argv[2])
    min_date = int(str(min_month) + "00")
    max_date = int(str(max_month) + "32")

    mode = sys.argv[3]
    if mode == 'rescor':
        program_pattern = "python get_rescor_by_shot_nomask.py {}\n"
        filename = "run_rescor"
    elif mode == "correction":
        filename = "run_correction"
        program_pattern = "python get_improved_ffskysub.py {} 1 0\n"

    survey = Survey(HDRVersion).return_astropy_table()
    here = (survey['date'] > min_date) * (survey['date'] < max_date)
    survey = survey[here]

    shotlist = survey['shotid']

    with open(filename, "w") as ff:
        for shotid in shotlist:
            ff.write(program_pattern.format(shotid))

    print(f"Finished writing to {filename}.")