from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from app.api.endpoints import getFingerPrint_bySMILES

app = FastAPI(
    title="Chemical Preprocessing",
    description="API for chemical SMILES standardization and fingerprint generation",
    version="0.1.0",
)

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(getFingerPrint_bySMILES.router, prefix="/getFingerPrint_bySMILES")

@app.get("/", tags=["Health"])
async def root():
    """Health check endpoint."""
    return {"status": "healthy", "message": "Chemical Preprocessing API is running"}


if __name__ == "__main__":
    import uvicorn

    uvicorn.run("app.main:app", host="0.0.0.0", port=8008, reload=True)
