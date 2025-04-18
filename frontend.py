import streamlit as st
import httpx
import os

API_URL = os.getenv("API_URL", "http://localhost:8008")
#API_URL = "http://localhost:8008/getFingerPrint_bySMILES"  

st.set_page_config(page_title="Chemical Fingerprint Generator", layout="centered")
st.title("SMILES to Fingerprint")

# Input field
smiles_input = st.text_input("Enter SMILES string:", placeholder="e.g., CC(=O)OC1=CC=CC=C1C(=O)O")

if st.button("Generate Fingerprint"):
    if not smiles_input.strip():
        st.warning("Please enter a valid SMILES string.")
    else:
        with st.spinner("Calling API..."):
            try:
                response = httpx.get(f"{API_URL}/getFingerPrint_bySMILES/{smiles_input}", timeout=10)

                if response.status_code == 200:
                    data = response.json()
                    st.success("Fingerprint generated successfully!")
                    st.text_area("Standardized SMILES", data["standardized_smiles"], height=100)

                    # Format fingerprint into rows of 64 for easier viewing
                    bitstring = data["fingerprint"]
                    bit_chunks = [bitstring[i:i+64] for i in range(0, len(bitstring), 64)]


                    st.markdown("### Fingerprint:")

                    styled_fingerprint = ""
                    for row in bit_chunks:
                        styled_row = "".join(
                            f"<span style='color:red;font-weight:bold;'>1</span>" if bit == 1
                            else "<span style='color:#90ee90;'>0</span>"
                            for bit in row
                        )
                        styled_fingerprint += styled_row + "<br>"

                    st.markdown(
                        f"""
                        <div style="width: 100%; border: 1px solid #ccc; border-radius: 6px;
                                    font-family: monospace; font-size: 1.2rem;
                                    background-color: #f9f9f9; overflow-x: auto;
                                    white-space: pre-wrap; padding: 0.25rem 0.25rem 0.25rem 0.25rem;
                                    line-height: 1.6rem; margin-bottom: 0;">
                            <pre style="margin: 0; padding: 0;">{styled_fingerprint}</pre>
                        </div>
                        """,
                        unsafe_allow_html=True,
                    )
                    

 
                elif response.status_code == 500:
                    error_msg = response.json().get("detail", "Unknown server error.")
                    st.error(f"Server Error: {error_msg}")
                else:
                    st.error(f"Unexpected Error: Status {response.status_code}")

            except Exception as e:
                st.error(f"Connection error: {e}")
