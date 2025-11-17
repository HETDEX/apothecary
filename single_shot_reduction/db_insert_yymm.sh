## example insert script
## change the paths and YYMM as needed
set -o noglob ; python3 make_report_db.py --db_name /scratch/03261/polonius/parallel/detect/image_db/elixer_reports_YYMM.db --img_dir /scratch/03261/polonius/parallel/all_pngs_YYYYMM --img_name YYMM*[0-9].png &
set -o noglob ; python3 make_report_db.py --db_name /scratch/03261/polonius/parallel/detect/image_db/elixer_reports_YYMM_nei.db --img_dir /scratch/03261/polonius/parallel/all_pngs_YYYYMM --img_name YYMM*_nei.png &
