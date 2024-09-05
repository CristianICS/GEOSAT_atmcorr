import os
import sys

def group_args(args: list):
    """Group command line params and its values inside a python dict."""
    args_dict = {}
    for arg in args:
        # Each argument should be '-arg=value'
        name, value = arg.split('=')
        args_dict[name] = value
    
    return args_dict

if __name__ == "__main__":
    # Check path for atmcorr utils
    ROOT = os.path.dirname(sys.argv[0])
    sys.path.append(ROOT)
    import atmcorr_utils as atm # type: ignore

    # Remove the script path and remain only the custom arguments
    sys.argv.pop(0)
    # Save a dictionary with the custom arguments
    args = group_args(sys.argv)

    # Open images
    try:
        main_folder = os.path.dirname(os.path.abspath(args['-folder']))
        img_folders = [args['-folder']]
    except:
        try:
            main_folder = os.path.abspath(args['-folders'])
            img_folders = os.listdir(main_folder)
        except:
            raise ValueError('The "-folder"/"-folders" param is missing or invalid.')
    
    # Perform the correction
    try:
        atm_key = args['-atm_key']
        geecloud_id = args['-geecloud_id'] if '-geecloud_id' in args.keys() else None
        try:
          for folder in img_folders:
            try:
              img = atm.Geosat(os.path.join(main_folder, folder))
              print(f'Compute image folder {os.path.basename(folder)}')
              img.atm_corr(atm_key, **{'gee_id':geecloud_id})
            except ValueError as ve:
              print(ve)
        except ValueError as ve:
            print(ve)
    except:
        raise ValueError(" ".join([
            'The "-atm_key" parameters is missing or invalid. Insert one of',
            'these: [ARC, DOS, COST, 6S]. Note that for DOS, COST and 6S',
            'corrections you must have an image with the ARC correction',
            'performed. The 6S correction require the "-geecloud_id" param.'
        ]))
