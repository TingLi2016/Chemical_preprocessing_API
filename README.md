# Chemical Preprocessing API

A containerized REST API for chemical SMILES, built with FastAPI. This API allows users to upload SMILES and perform chemical cleaning and fingerprint generation. 


## Features

- **Fingerprint Generation**: Support chemical cleaning and fingerprint generation

## Installation

### Prerequisites

- Docker and Docker Compose (recommended)
- Python 3.10+ (for local development)
- uv package manager (for local development)

### Option 1: Using Docker (Recommended)

1. Clone the repository:
   ```bash
   git clone <repository-url>
   cd chemical_preprocessing
   ```

2. Build and run the containers:
   ```bash
   docker-compose up -d
   ```

3. Access the application:
   - Frontend: http://localhost:8502
   - Backend API: http://localhost:8008
   - API Documentation (Swagger UI): http://localhost:8008/docs
   - API Documentation (ReDoc): http://localhost:8008/redoc

### Option 2: Local Development

1. Clone the repository:
   ```bash
   git clone <repository-url>
   cd chem-preprocessing
   ```

2. Install dependencies using uv:
   ```bash
   uv pip install -e ".[dev]"
   ```

3. Run the backend:
   ```bash
   uv run uvicorn app.main:app --host 0.0.0.0 --port 8000 --reload
   ```

4. Run the frontend (in a separate terminal):
   ```bash
   uv run streamlit run frontend.py -- --server.port=8501 --server.address=0.0.0.0
   ```

5. Access the application:
   - Frontend: http://localhost:8502
   - Backend API: http://localhost:8008
   - API Documentation (Swagger UI): http://localhost:8008/docs
   - API Documentation (ReDoc): http://localhost:8008/redoc

## Usage

### Using the Dashboard

1. **getFingerPrint_bySMILES**:
   - Input SMILES string

## API Endpoints

### API Documentation
- Interactive API docs (Swagger UI): http://localhost:8008/docs
- Alternative API docs (ReDoc): http://localhost:8008/redoc

### getFingerPrint_bySMILES Endpoints
- `GET /`: Get the SMILES from user input


## Development

### Running Tests

Using Docker:
```bash
docker build --target test .
```

Using uv:
```bash
uv run pytest -v
```

### Project Structure

```
chemical_preprocessing_api/
├── app/                      # Backend API
│   ├── api/                  # API endpoints
│   │   └── endpoints/        # Endpoint modules
│   ├── models/               # Data models
│   ├── services/             # Analysis services
├── tests/                    # Test suite
├── frontend.py               # Streamlit frontend
├── Dockerfile                # Multi-stage Docker build
├── docker-compose.yml        # Docker Compose configuration
└── environment.yml            # Conda environment file, made this special for RDKit package
└── pyproject.toml            # Project dependencies
```

## Technical Details

### Backend (FastAPI)
- Built with FastAPI for high-performance API endpoints
- Uses Pydantic for data validation
- Implements chemical cleaning and fingerprint generation for SMILES

### Frontend (Streamlit)
- Built with Streamlit for rapid dashboard development


## License

[MIT License](LICENSE)
