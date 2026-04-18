import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

export default defineConfig({
  plugins: [react()],
  base: '/brcapath-rx/confounding-demo/',
  build: {
    outDir: '../../docs/confounding-demo',
    emptyOutDir: true,
  },
})
