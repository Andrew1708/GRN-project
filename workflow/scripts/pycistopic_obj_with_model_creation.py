import os
import pickle
from pycisTopic.lda_models import run_cgs_models_mallet
from pycisTopic.lda_models import evaluate_models
import argparse

n_topics = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
os.environ['MALLET_MEMORY'] = '200G'

def create_cistopic_object_with_model(pkl_file_path, mallet_path, temp):
    # Load the cistopic object
    with open(pkl_file_path,'rb') as f:
        cistopic_object = pickle.load(f)

    # Run MALLET model
    models = run_cgs_models_mallet(
        cistopic_object,
        n_topics=n_topics,
        n_cpu=len(n_topics),
        random_state=555,
        alpha=50,
        alpha_by_topic=True,
        eta=0.1,
        eta_by_topic=False,
        tmp_path= temp,
        mallet_path=mallet_path,
    )

    # Select the best model
    model = evaluate_models(
        models,
        return_model = True,
    )

    # Add model to cistopic object
    cistopic_object.add_LDA_model(model)

    return cistopic_object

def parser():
    parser = argparse.ArgumentParser(description="Process directories for barcode and region files.")
    parser.add_argument("--cistopic_path", required=True, help="Path to the cistopic object .pkl file")
    parser.add_argument("--temp_dir", required=True, help="Temporary directory for models")
    parser.add_argument("--out_dir", required=True, help="Output directory for models")
    parser.add_argument("--mallet_path", required=True, help="Path to the Mallet binary")
    return parser.parse_args()


def main():
    args = parser()

    cistopic_path = args.cistopic_path
    temp = args.temp_dir
    out_dir = args.out_dir
    mallet_path = args.mallet_path

    # (Optional) Create directories if needed
    os.makedirs(temp, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)

    # Create cistopic object with model
    cistopic_object = create_cistopic_object_with_model(cistopic_path, mallet_path, temp)
    # Save the updated cistopic object
    out_cistopic_file = os.path.join(out_dir, f"{cistopic_object.project}_model_pycistopic.pkl")
    with open(out_cistopic_file, "wb") as f:
        pickle.dump(cistopic_object, f)

    print(f"Saved: {out_cistopic_file}")


if __name__ == "__main__":
    main()